#pragma once

#include "measurement.h"
#include "state.h"
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <numeric>
#include <pb_utils/numbers.h>
#include <pb_utils/timer.h>
#include <ranges>

namespace {
struct LinearSystem {
  /* Eigen::MatrixXd H; */
  std::vector<Eigen::Triplet<double>> jacobian_triplets;
  std::vector<Eigen::Triplet<double>> hessian_triplets;
  Eigen::SparseMatrix<double> J;
  Eigen::SparseMatrix<double> H;
  Eigen::VectorXd g;
  double chi_square{0.0};
  bool sp_pattern_analyzed{false};
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sp_solver;
  /* Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> sp_solver; */

  LinearSystem(const long system_size, const long measurement_size) : g(Eigen::VectorXd::Zero(system_size)) {
    /* allocates memory for jacobian and hessian triples, the sparse jacobian and
     * hessian matrices */
    jacobian_triplets.reserve(measurement_size * 2);
    hessian_triplets.reserve(measurement_size * 4);
    Eigen::SparseMatrix<double> J(4 * measurement_size, system_size);
    Eigen::SparseMatrix<double> H(system_size, system_size);
  }

  /* clears all the entries to allow reusing the allocated memory */
  void clear() {
    jacobian_triplets.clear();
    hessian_triplets.clear();
    J.setZero();
    H.setZero();
    g.setZero();
    chi_square = 0.0;
    sp_pattern_analyzed = false;
  }
  void convert_triplets() {
    J.setFromTriplets(jacobian_triplets.begin(), jacobian_triplets.end());
    H.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end());
  }
  /* calculate the gauss newton step */
  auto solve() {
    const auto system_size = H.outerSize();
    Eigen::VectorXd dx = Eigen::VectorXd::Zero(system_size);
    const auto& H_block = H.bottomRightCorner(system_size - 1, system_size - 1);
    const auto& g_block = g.tail(system_size - 1);
    if (!sp_pattern_analyzed) {
      /* split the compute step since our sparsity structure remains the same accross iterations */
      sp_solver.analyzePattern(H_block);
      sp_pattern_analyzed = true;
    }
    sp_solver.factorize(H_block);
    dx.tail(system_size - 1) = sp_solver.solve(g_block);
    return dx;
  }
};
} // namespace

namespace pidgeot {
/* powell's */
/* see Methods for Non-Linear Least Squares Problems, Madsen et al. 2004 */
class DogLegSolver {
private:
  int _max_iter;
  State _state;
  Measurement _measurement;
  double _chi_square_thresh;
  double _dx_sqnorm_thresh;

public:
  DogLegSolver(int max_iter,
               const State& initial_state,
               const Measurement& measurement,
               double chi_square_thresh = 1e-8,
               double dx_sqnorm_thresh = 1e-16)
      : _state(initial_state),
        _measurement(measurement),
        _max_iter(max_iter),
        _chi_square_thresh(chi_square_thresh),
        _dx_sqnorm_thresh(dx_sqnorm_thresh) {}

  auto solve(bool verbose = false) {
    const long system_size = pb_utils::saturate_cast<long>(_state.size());
    const long measurement_size = pb_utils::saturate_cast<long>(_measurement.size());

    pb_utils::Timer lsq_timer("LSQ Optimization");

    LinearSystem linear_system(system_size, measurement_size);
    /* Eigen::VectorXd gn_step = Eigen::VectorXd::Zero(system_size - 1); */
    /* Eigen::VectorXd sd_step = Eigen::VectorXd::Zero(system_size - 1); */
    /* Eigen::VectorXd dl_step = Eigen::VectorXd::Zero(system_size - 1); */

    // lambda
    auto rotation_derived = [](const double angle) {
      Eigen::Matrix2d R_dot{
          {-std::sin(angle), -std::cos(angle)},
          { std::cos(angle), -std::sin(angle)}
      };
      return R_dot;
    };
    // need to catch some stuff by reference before this lambda
    auto update_linear_system = [&](const auto& enumerated_measurement) {
      const auto& [z_idx, z] = enumerated_measurement;
      const auto from_idx = z.from_state_idx;
      const auto to_idx = z.to_state_idx;
      const auto& from = _state[from_idx].rotation;
      const auto& to = _state[to_idx].rotation;
      Eigen::Vector4d h_ij = (from.inverse() * to).toRotationMatrix().reshaped();
      Eigen::Vector4d e = h_ij - z.rotation.toRotationMatrix().reshaped();
      Eigen::Vector4d Ji = (rotation_derived(from.angle()).transpose() * to.toRotationMatrix()).reshaped();
      Eigen::Vector4d Jj = (from.inverse().toRotationMatrix() * rotation_derived(to.angle())).reshaped();
      const auto idx_range = std::views::iota(z_idx, z_idx + 4);
      const auto zipped_Ji_entries = std::views::zip(idx_range, std::views::repeat(from_idx, 4), Ji);
      const auto zipped_Jj_entries = std::views::zip(idx_range, std::views::repeat(to_idx, 4), Jj);
      linear_system.jacobian_triplets.append_range(zipped_Ji_entries);
      linear_system.jacobian_triplets.append_range(zipped_Jj_entries);
      linear_system.hessian_triplets.emplace_back(from_idx, from_idx, Ji.transpose() * Ji);
      linear_system.hessian_triplets.emplace_back(to_idx, to_idx, Jj.transpose() * Jj);
      linear_system.hessian_triplets.emplace_back(from_idx, to_idx, Ji.transpose() * Jj);
      linear_system.hessian_triplets.emplace_back(to_idx, from_idx, Ji.transpose() * Jj);
      linear_system.g(from_idx) += -Ji.transpose() * e;
      linear_system.g(to_idx) += -Jj.transpose() * e;
      linear_system.chi_square += e.squaredNorm();
    };

    int iter = 0;
    while (iter < _max_iter) {
      lsq_timer.tick();
      linear_system.clear();
      std::ranges::for_each(std::views::enumerate(_measurement), update_linear_system);
      linear_system.convert_triplets();
      auto dx = linear_system.solve();
      _state.box_plus(dx);
      auto dx_sqnorm = dx.squaredNorm();

      auto chi_square = linear_system.chi_square;
      std::cout << "Iter: " << iter << "/" << _max_iter - 1 << " and chi_squared = " << chi_square
                << " and delta_x (sq. norm) = " << dx_sqnorm << ", took " << lsq_timer.tock() << "s\n";
      if (verbose) {
        std::cout << "g is\n" << linear_system.g.transpose() << " \n";
      }
      if (chi_square < _chi_square_thresh || dx_sqnorm < _dx_sqnorm_thresh) {
        break;
      }
      iter++;
    }
    return _state;
  }
};
}; // namespace pidgeot
