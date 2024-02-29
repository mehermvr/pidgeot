#include "dogLegSolver.h"
#include "steepestDescentSolver.h"
#include <cmath>
#include <stdexcept>

namespace pidgeot {

State DogLegSolver::solve(bool verbose) {
  pb_utils::Timer lsq_timer("Dog Leg Optimization");

  const long system_size = pb_utils::saturate_cast<long>(_x.size());
  const long measurement_size = pb_utils::saturate_cast<long>(_measurement.size());

  LinearSystem linear_system(system_size, measurement_size);

  int iter = 0;
  while (iter < _max_iter) {
    if (verbose) {
      lsq_timer.tick();
    }
    linear_system.build(_x, _measurement);
    const Eigen::VectorXd dx_gn = linear_system.solve();
    if (linear_system.chi_square <= _eps_3) {
      /* termination 3 satisfied */
      /* state is not updated with dx_dl */
      if (verbose) {
        std::cout << "chi_square is less than threshold. terminating\n";
      }
      break;
    }
    // to use in logging later
    const double g_infinite_norm = linear_system.g.lpNorm<Eigen::Infinity>();
    if (g_infinite_norm <= _eps_1) {
      /* termination 1 satisfied */
      /* state is not updated with dx_dl */
      if (verbose) {
        std::cout << "max elem in the gradient is less than threshold. terminating\n";
      }
      break;
    }
    double alpha{}; // modified in solve() call
    const Eigen::VectorXd dx_sd = SteepestDescentSolver::solve(linear_system, alpha);
    const auto& [dx_dl, dl_case, beta] = dogleg_step(dx_sd, dx_gn);
    const auto dx_dl_norm = dx_dl.norm();
    if (dx_dl_norm <= _eps_2) {
      /* termination 2 satisfied */
      /* state is not updated with dx_dl */
      if (verbose) {
        std::cout << "dx_dl dogleg_step norm is less than threshold. terminating\n";
      }
      break;
    }
    State x_new = _x;
    x_new.box_plus(dx_dl);
    const double gain_ratio = calculate_gain_ratio(x_new, linear_system, dl_case, alpha, beta);
    if (gain_ratio > 0) {
      _x = x_new;
      /* leave the termination criteria to the next iteration when the system is solved */
    }
    if (gain_ratio > 0.75) {
      _trust_radius = std::max(_trust_radius, 3 * dx_dl_norm);
    } else if (gain_ratio < 0.25) {
      _trust_radius = _trust_radius / 2;
    }
    // else trust radius is unchanged

    if (verbose) {
      /* std::cout << "Gain ratio: " << gain_ratio << " and trust radius: " << _trust_radius << "\n"; */
      std::cout << "Iter: " << iter << "/" << _max_iter - 1 << " - T1: ||g||_inf = " << g_infinite_norm
                << ", T2: ||dx||_abs_change = " << dx_dl_norm << ", T3: chi_square = " << linear_system.chi_square
                << ", took " << lsq_timer.tock() << "s\n";
    }
    iter++;
  }
  return _x;
}
std::tuple<Eigen::VectorXd, int, double> DogLegSolver::dogleg_step(const auto& sd_step, const auto& gn_step) {
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
double DogLegSolver::calculate_gain_ratio(const State& x_new,
                                          const LinearSystem& linear_system,
                                          const int dl_case,
                                          const double alpha,
                                          const double beta) const {
  double linear_model_gain{};
  const double F = 0.5 * linear_system.chi_square;
  switch (dl_case) {
  case 1:
    /* std::cout << "dl case 1: GN step\n"; */
    // F(x), where F(x) = 0.5 f(x).T * f(x). for us f(x) is e(x). and e.T * e is stored in chi_square
    linear_model_gain = F;
    break;
  case 2:
    /* std::cout << "dl case 2: trust radius length in gradient direction\n"; */
    // alpha * linear_system.g is just sd_step. can capture that and reduce some flops
    linear_model_gain = _trust_radius * (2 * (alpha * linear_system.g).norm() - _trust_radius) / (2 * alpha);
    break;
  case 3:
    /* std::cout << "dl case 3: trust radius length in between gradient and gauss newton direction\n"; */
    linear_model_gain = 0.5 * alpha * (1 - beta) * (1 - beta) * linear_system.g.squaredNorm() + beta * (2 - beta) * F;
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
