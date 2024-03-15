#include "dogLegSolver.h"
#include "linearSystem.h"
#include <cmath>
#include <stdexcept>
#include <utility>

using DiagonalMatrixXd = Eigen::DiagonalMatrix<double, Eigen::Dynamic>;
namespace pidgeot {

void DogLegSolver::_calculate_trust_diagonal(const LinearSystem& linear_system) {
  _trust_diagonal = Eigen::VectorXd::Ones(linear_system.g.size());
  /* _trust_diagonal = linear_system.H.diagonal(); */
  /* TODO: test the effect */
  _trust_diagonal = _trust_diagonal.cwiseSqrt();
  // these are values from ceres defaults
  // https://github.com/ceres-solver/ceres-solver/blob/62c03d6ff3b1735d4b0a88006486cfde39f10063/docs/source/nnls_solving.rst#L1517
  double min_diag_val = 1e-6;
  double max_diag_val = 1e32;
  for (auto& diag_val : _trust_diagonal) {
    diag_val = std::min(max_diag_val, std::max(min_diag_val, diag_val));
  }
}
Eigen::VectorXd DogLegSolver::_trust_scale_vector(const Eigen::VectorXd& vec) const {
  Eigen::VectorXd scaled_vec = _trust_diagonal.array() * vec.array();
  return scaled_vec;
}
Eigen::VectorXd DogLegSolver::_trust_unscale_vector(const Eigen::VectorXd& scaled_vec) const {
  Eigen::VectorXd vec = _trust_diagonal.array().inverse() * scaled_vec.array();
  return vec;
}
Eigen::VectorXd DogLegSolver::_calculate_gauss_newton_step(const LinearSystem& linear_system) {
  const auto& H = linear_system.H;
  const auto& g = linear_system.g;
  /* the system is (H:= J.T*J) * dx = -g */
  const auto system_size = H.outerSize();
  if (!_sparse_pattern_analyzed) {
    /* split the compute step since our sparsity structure remains the same accross iterations */
    _sparse_solver.analyzePattern(H);
    _sparse_pattern_analyzed = true;
  }
  _sparse_solver.factorize(H);
  Eigen::VectorXd dx = Eigen::VectorXd::Zero(system_size);
  dx = _sparse_solver.solve(-1 * g);
  // scaling the system essentially leads
  // D.inv() * H * D.inv() * dx_scaled = - D.inv() * g
  // H * D.inv() * dx_scaled = - g
  // D.inv() * dx_scaled = dx_gn
  // dx_scaled = D * dx_gn, which makes sense
  _trust_scale_vector(dx);
  return dx;
}
std::pair<double, Eigen::VectorXd> DogLegSolver::_calculate_cauchy_step(const LinearSystem& linear_system) {
  // the scaled system, H is D.inv() * H * D.inv() and g is D.inv() * g.
  // calculating alpha involves norm of scaled_g and scaled_g.T * scaled_H * scaled_g
  // so g.T * D.inv().T * D.inv() * H * D.inv() * D.inv() * g
  // => (D.inv() * D.inv() * g).T * H * (D.inv() * D.inv() * g)
  // so instead of scaling H, you can get away with double scaling g
  const auto& g = linear_system.g;
  Eigen::VectorXd scaled_g = _trust_unscale_vector(g);
  Eigen::VectorXd double_scaled_g = _trust_unscale_vector(scaled_g);
  const auto& H = linear_system.H;
  const double alpha = scaled_g.squaredNorm() / (double_scaled_g.transpose() * H * double_scaled_g);
  Eigen::VectorXd sd_step = -alpha * scaled_g;
  return std::make_pair(alpha, sd_step);
}

State DogLegSolver::solve(bool verbose) {
  pb_utils::Timer lsq_timer("Dog Leg Optimization");

  auto termination_criteria_1 = [&](const double g_infinite_norm) {
    if (g_infinite_norm <= _eps_1) {
      if (verbose) {
        std::cout << "max elem in the gradient is less than threshold. terminating\n";
      }
      return true;
    }
    return false;
  };
  auto termination_criteria_2 = [&](const double dx_dl_norm) {
    if (dx_dl_norm <= _eps_2) {
      if (verbose) {
        std::cout << "dx_dl dogleg_step norm is less than threshold. terminating\n";
      }
      return true;
    }
    return false;
  };
  auto termination_criteria_3 = [&](const double chi_square) {
    if (chi_square <= _eps_3) {
      if (verbose) {
        std::cout << "chi_square is less than threshold. terminating\n";
      }
      return true;
    }
    return false;
  };

  const long system_size = pb_utils::saturate_cast<long>(_x.size());
  const long measurement_size = pb_utils::saturate_cast<long>(_measurement.size());

  LinearSystem linear_system(system_size, measurement_size);
  int iter = 0;
  while (iter < _max_iter) {
    if (verbose) {
      lsq_timer.tick();
    }
    linear_system.build(_x, _measurement);
    // to use in logging later
    const double g_infinite_norm = linear_system.g.lpNorm<Eigen::Infinity>();

    _calculate_trust_diagonal(linear_system);
    // scaled by trust region
    const Eigen::VectorXd dx_gn = _calculate_gauss_newton_step(linear_system);
    const auto& [alpha, dx_sd] = _calculate_cauchy_step(linear_system);
    const auto& [sdx_dl, dl_case, beta] = _calculate_dogleg_step(dx_sd, dx_gn);
    Eigen::VectorXd dx_dl = _trust_unscale_vector(sdx_dl);
    const auto dx_dl_norm = dx_dl.norm();
    if (termination_criteria_1(g_infinite_norm) || termination_criteria_2(dx_dl_norm) ||
        termination_criteria_3(linear_system.chi_square)) {
      break;
    }
    State x_new = _x;
    x_new.box_plus(dx_dl);
    const double gain_ratio = _calculate_gain_ratio(x_new, linear_system, dl_case, alpha, beta);

    if (gain_ratio > 0) {
      _x = x_new;
    }
    if (gain_ratio > 0.75) {
      _trust_radius = std::max(_trust_radius, 3 * dx_dl_norm);
    } else if (gain_ratio < 0.25) {
      _trust_radius = _trust_radius / 2;
    }
    // else trust radius is unchanged

    if (verbose) {
      /* std::cout << "Gain ratio: " << gain_ratio << " and trust radius: " << _trust_radius << "\n"; */
      std::cout << "Iter: " << iter << "/" << _max_iter - 1 << "(DL: " << dl_case
                << ") - T1: ||g||_inf = " << g_infinite_norm << ", T2: ||dx||_abs_change = " << dx_dl_norm
                << ", T3: chi_square = " << linear_system.chi_square << ", radius: " << _trust_radius << ", took "
                << lsq_timer.tock() << "s\n";
    }
    iter++;
  }
  return _x;
}
std::tuple<Eigen::VectorXd, int, double> DogLegSolver::_calculate_dogleg_step(const auto& sd_step,
                                                                              const auto& gn_step) {
  const auto trust_radius_sq = _trust_radius * _trust_radius;
  /* this is the step size along gn - sd direction */
  auto calculate_beta = [&]() {
    /* renaming to make life easy to follow madsen, page 31 */
    const auto& a = sd_step;
    const auto& b = gn_step;
    const double c = a.transpose() * (b - a);
    double beta = 0.0;
    const double a_sqnorm = a.squaredNorm();
    const double b_minus_a_sqnorm = (b - a).squaredNorm();
    const double sqrt_term = std::sqrt(c * c + b_minus_a_sqnorm * (trust_radius_sq - a_sqnorm));
    if (c <= 0) {
      beta = (-c + sqrt_term) / b_minus_a_sqnorm;
    } else {
      beta = (trust_radius_sq - a_sqnorm) / (c + sqrt_term);
    }
    return beta;
  };

  double beta = -1.0;
  if (gn_step.squaredNorm() <= trust_radius_sq) {
    return std::make_tuple(gn_step, 1, beta);
  } else if (sd_step.squaredNorm() >= trust_radius_sq) {
    Eigen::VectorXd dl_step = _trust_radius * sd_step / sd_step.norm();
    return std::make_tuple(dl_step, 2, beta);
  } else {
    beta = calculate_beta();
    Eigen::VectorXd dl_step = sd_step + beta * (gn_step - sd_step);
    return std::make_tuple(dl_step, 3, beta);
  }
};
double DogLegSolver::_calculate_gain_ratio(const State& x_new,
                                           const LinearSystem& linear_system,
                                           const int dl_case,
                                           const double alpha,
                                           const double beta) const {
  double linear_model_gain{};
  const double F = 0.5 * linear_system.chi_square;
  Eigen::VectorXd scaled_g = _trust_unscale_vector(linear_system.g);
  switch (dl_case) {
  case 1:
    /* std::cout << "dl case 1: GN step\n"; */
    // F(x), where F(x) = 0.5 f(x).T * f(x). for us f(x) is e(x). and e.T * e is stored in chi_square
    linear_model_gain = F;
    break;
  case 2:
    /* std::cout << "dl case 2: trust radius length in gradient direction\n"; */
    // alpha * linear_system.g is just sd_step. can capture that and reduce some flops
    linear_model_gain = _trust_radius * (2 * (alpha * scaled_g).norm() - _trust_radius) / (2 * alpha);
    break;
  case 3:
    /* std::cout << "dl case 3: trust radius length in between gradient and gauss newton direction\n"; */
    linear_model_gain = 0.5 * alpha * (1 - beta) * (1 - beta) * scaled_g.squaredNorm() + beta * (2 - beta) * F;
    break;
  default:
    throw std::invalid_argument("should never reach");
  }
  // linear model gain is calculated. need non linear model gain now
  auto calculate_F_new = [&]() {
    double error{};
    std::ranges::for_each(_measurement, [&](const AtomicMeasurement& z) {
      /* const auto& [from, to] = z.get_associated_states(x_new); */
      const auto& from = x_new[z.from_state_idx];
      const auto& to = x_new[z.to_state_idx];
      auto e = z.calculate_error(from, to);
      error += e.squaredNorm();
    });
    return 0.5 * error;
  };
  const double F_new = calculate_F_new();
  const double non_linear_model_gain = F - F_new;
  return non_linear_model_gain / linear_model_gain;
}
} // namespace pidgeot
