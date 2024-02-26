#pragma once

#include "measurement.h"
#include "state.h"
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <cmath>
#include <pb_utils/numbers.h>
#include <pb_utils/timer.h>
#include <stdexcept>

namespace pidgeot {
struct LinearSystem {
  /* Eigen::MatrixXd H; */
  std::vector<Eigen::Triplet<double>> hessian_triplets;
  /* hessian, but more like J.T * J, ignores the second order terms. if the error at optimum is 0, then H = J.T * J and
   * you get quadratic convergence with GN. see madsen */
  Eigen::SparseMatrix<double> H;
  /* the gradient. the system is (J.T*J)h = -g */
  Eigen::VectorXd g;
  double chi_square{0.0};
  bool sp_pattern_analyzed{false};
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sp_solver;
  /* Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> sp_solver; */

  LinearSystem(const long system_size, const long measurement_size)
      : g(Eigen::VectorXd::Zero(system_size)), H(system_size, system_size) {
    /* allocates memory for hessian triples, the sparse hessian matrices */
    hessian_triplets.reserve(measurement_size * 4);
  }

  /* clears all the entries to allow reusing the allocated memory */
  void clear() {
    hessian_triplets.clear();
    H.setZero();
    g.setZero();
    chi_square = 0.0;
    sp_pattern_analyzed = false;
  }
  void convert_triplets() { H.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end()); }

  /* calculate the gauss newton step */
  auto solve() {
    /* the system is (H:= J.T*J) * dx = -g */
    const auto system_size = H.outerSize();
    Eigen::VectorXd dx = Eigen::VectorXd::Zero(system_size);
    /* TODO: this is because we fix the first state. there has to be a better way */
    const auto& H_block = H.bottomRightCorner(system_size - 1, system_size - 1);
    const auto& g_block = g.tail(system_size - 1);
    if (!sp_pattern_analyzed) {
      /* split the compute step since our sparsity structure remains the same accross iterations */
      sp_solver.analyzePattern(H_block);
      sp_pattern_analyzed = true;
    }
    sp_solver.factorize(H_block);
    dx.tail(system_size - 1) = sp_solver.solve(-g_block);
    return dx;
  }
};

/* see Methods for Non-Linear Least Squares Problems, Madsen et al. 2004 */
/* Termination Criteria - Powell's Dog Leg Solver
 * 1. inifinity_norm(g) <= eps_1 => gradient is close to 0, we are at a local minimum
 * 2. (dx).norm() <= eps_2 => absolute change in x is small. NOTE this is different from madsen's termination criterion
 * which is based on relative change
 * 3. iter >= iter_max => prevent inifinite iterations
 * 4. inifinity_norm(error) <= eps_3 => error is close to 0. We do a threshold on chi_square instead since this is
 * easier in the current implementation
 * 5. trust_radius <= eps_2(x.norm + eps_2) => because if this is true, condition 2 is gauranteed. NOTE this is not
 * implemented since this also involves x.norm like criterion 2 which is a pain. its gauranteed, let the next iteration
 * handle it
 */
class DogLegSolver {
private:
  int _max_iter;
  State _x;
  Measurement _measurement;
  double _eps_1;
  double _eps_2;
  double _eps_3;
  double _trust_radius;

public:
  /*
   * param eps_1: if max element in gradient at an iteration is less than eps_1, terminate
   * param eps_2: if obsolute change in x across an interation (dx) is less than eps_2, terminate. this is different
   * from madsen, this the absolute is easier to implement than relative. (state_norm is a pain)
   * param eps_3: if chi_square in an iteration is less than eps_3, terminate. this is different from Madsen
   */
  DogLegSolver(int max_iter,
               const State& initial_state,
               const Measurement& measurement,
               double trust_radius = 1,
               double eps_1 = 1e-16,
               double eps_2 = 1e-20,
               double eps_3 = 1e-16)
      : _x(initial_state),
        _measurement(measurement),
        _max_iter(max_iter),
        _trust_radius(trust_radius),
        _eps_1(eps_1),
        _eps_2(eps_2),
        _eps_3(eps_3) {}

  auto solve(bool verbose = false) {
    pb_utils::Timer lsq_timer("LSQ Optimization");

    const long system_size = pb_utils::saturate_cast<long>(_x.size());
    const long measurement_size = pb_utils::saturate_cast<long>(_measurement.size());

    LinearSystem linear_system(system_size, measurement_size);

    // its lambda-ing time
    auto rotation_derived = [](const double angle) {
      Eigen::Matrix2d R_dot{
          {-std::sin(angle), -std::cos(angle)},
          { std::cos(angle), -std::sin(angle)}
      };
      return R_dot;
    };
    /* this is reused later in calculating the non-linear model gain */
    auto calculate_error_for_one_z = [](const State& x, const AtomicMeasurement& z) {
      const auto from_idx = z.from_state_idx;
      const auto to_idx = z.to_state_idx;
      const auto& from = x[from_idx].rotation;
      const auto& to = x[to_idx].rotation;
      Eigen::Vector4d h_ij = (from.inverse() * to).toRotationMatrix().reshaped();
      Eigen::Vector4d e = h_ij - z.rotation.toRotationMatrix().reshaped();
      return std::make_tuple(e, from_idx, to_idx, from, to);
    };
    // need to catch some stuff by reference before this lambda
    auto update_linear_system = [&](const auto& z) {
      const auto& [e, from_idx, to_idx, from, to] = calculate_error_for_one_z(_x, z);
      Eigen::Vector4d Ji = (rotation_derived(from.angle()).transpose() * to.toRotationMatrix()).reshaped();
      Eigen::Vector4d Jj = (from.inverse().toRotationMatrix() * rotation_derived(to.angle())).reshaped();
      linear_system.hessian_triplets.emplace_back(from_idx, from_idx, Ji.transpose() * Ji);
      linear_system.hessian_triplets.emplace_back(to_idx, to_idx, Jj.transpose() * Jj);
      linear_system.hessian_triplets.emplace_back(from_idx, to_idx, Ji.transpose() * Jj);
      linear_system.hessian_triplets.emplace_back(to_idx, from_idx, Ji.transpose() * Jj);
      /* g is the gradient, not the rhs in the normal equations */
      linear_system.g(from_idx) += Ji.transpose() * e;
      linear_system.g(to_idx) += Jj.transpose() * e;
      linear_system.chi_square += e.squaredNorm();
    };
    auto rebuild_linear_system_mats = [&]() {
      linear_system.clear();
      std::ranges::for_each(_measurement, update_linear_system);
      linear_system.convert_triplets();
    };
    auto gauss_newton_step = [&]() { return linear_system.solve(); };
    /* needed when calculated steepest descent step. this is the step size */
    auto calculate_alpha_sd = [&]() {
      const auto& g = linear_system.g;
      const auto& H = linear_system.H;
      double alpha = g.squaredNorm() / (g.transpose() * H * g);
      return alpha;
    };
    auto steepest_descent_step = [&]() {
      const auto alpha = calculate_alpha_sd();
      /* first element will be 0, because in solving the linear system we also fix the first elem of g */
      Eigen::VectorXd sd_step = -alpha * linear_system.g;
      return std::make_tuple(sd_step, alpha);
    };
    /* needed when calculating dogleg step. this is the step size along gn - sd direction */
    auto calculate_beta_dl = [&](const auto& sd_step, const auto& gn_step) {
      const auto trust_radius_sq = _trust_radius * _trust_radius;
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
    /* returns dl_step, the dl switch case, and the beta value. beta is only meaningful in case 3 */
    /* we keep track of which dl case was the iteration, so we can calculate the linear model diff L(0) - L(dl_dx)
     * without keeping track of the Jacobian */
    auto dogleg_step = [&](const auto& sd_step, const auto& gn_step) {
      const auto trust_radius_sq = _trust_radius * _trust_radius;
      double beta = -1.0;
      if (gn_step.squaredNorm() <= trust_radius_sq) {
        return std::make_tuple(gn_step, 1, beta);
      } else if (sd_step.squaredNorm() >= trust_radius_sq) {
        Eigen::VectorXd dl_step = _trust_radius * sd_step / sd_step.norm();
        return std::make_tuple(dl_step, 2, beta);
      } else {
        beta = calculate_beta_dl(sd_step, gn_step);
        Eigen::VectorXd dl_step = sd_step + beta * (gn_step - sd_step);
        return std::make_tuple(dl_step, 3, beta);
      }
    };
    auto calculate_F = [&](const State& x) {
      double error{};
      std::ranges::for_each(_measurement, [&](const auto& z) {
        Eigen::Vector4d e;
        std::tie(e, std::ignore, std::ignore, std::ignore, std::ignore) = calculate_error_for_one_z(x, z);
        error += e.squaredNorm();
      });
      return 0.5 * error;
    };
    auto calculate_gain_ratio = [&](const State& x_new, const int dl_case, const double alpha, const double beta) {
      double linear_model_gain{};
      const double F = 0.5 * linear_system.chi_square;
      /* std::cout << "current F is " << F << "\n"; */
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
        linear_model_gain =
            0.5 * alpha * (1 - beta) * (1 - beta) * linear_system.g.squaredNorm() + beta * (2 - beta) * F;
        break;
      default:
        throw std::invalid_argument("should never reach");
      }
      // linear model gain is calculated
      const double F_new = calculate_F(x_new);
      /* std::cout << "new F is " << F_new << "\n"; */
      const double non_linear_model_gain = F - F_new;
      return non_linear_model_gain / linear_model_gain;
    };
    int iter = 0;
    while (iter < _max_iter) {
      if (verbose) {
        lsq_timer.tick();
      }
      rebuild_linear_system_mats();
      const Eigen::VectorXd dx_gn = gauss_newton_step();
      if (linear_system.chi_square <= _eps_3) {
        /* termination 3 satisfied */
        /* state is not updated with dx_dl */
        if (verbose) {
          std::cout << "chi_square is less than threshold. terminating";
        }
        break;
      }
      if (linear_system.g.maxCoeff() <= _eps_1) {
        /* termination 1 satisfied */
        /* state is not updated with dx_dl */
        if (verbose) {
          std::cout << "max elem in the gradient is less than threshold. terminating";
        }
        break;
      }
      const auto& [dx_sd, alpha] = steepest_descent_step();
      const auto& [dx_dl, dl_case, beta] = dogleg_step(dx_sd, dx_gn);
      /* auto dx_dl = dx_sd; */
      const auto dx_dl_norm = dx_dl.norm();
      if (dx_dl_norm <= _eps_2) {
        /* termination 2 satisfied */
        /* state is not updated with dx_dl */
        if (verbose) {
          std::cout << "dx_dl dogleg_step norm is less than threshold. terminating";
        }
        break;
      }
      State x_new = _x;
      x_new.box_plus(dx_dl);
      const double gain_ratio = calculate_gain_ratio(x_new, dl_case, alpha, beta);
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
        std::cout << "Gain ratio: " << gain_ratio << " and trust radius: " << _trust_radius << "\n";
        std::cout << "Iter: " << iter << "/" << _max_iter - 1 << " and chi_squared = " << linear_system.chi_square
                  << " and delta_x_dl (norm) = " << dx_dl_norm << ", took " << lsq_timer.tock() << "s\n";
      }
      iter++;
    }
    return _x;
  }
};
}; // namespace pidgeot
