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

namespace {
struct LinearSystemEntry {
  /* Eigen::MatrixXd H; */
  Eigen::VectorXd g;
  double chi_square{0.0};

  explicit LinearSystemEntry(const long size) : g(Eigen::VectorXd::Zero(size)) {}

  LinearSystemEntry& operator+=(const LinearSystemEntry& other) {
    this->g += other.g;
    this->chi_square += other.chi_square;
    return *this;
  }

  friend LinearSystemEntry operator+(LinearSystemEntry lhs, const LinearSystemEntry& rhs) {
    lhs += rhs;
    return lhs;
  }
};
} // namespace

namespace pidgeot {
class Solver {
private:
  int _max_iter;
  State _state;
  Measurement _measurement;
  double _chi_square_thresh;
  double _dx_sqnorm_thresh;

public:
  Solver(int max_iter,
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

    // lambda
    auto rotation_derived = [](const double angle) {
      Eigen::Matrix2d R_dot{
          {-std::sin(angle), -std::cos(angle)},
          { std::cos(angle), -std::sin(angle)}
      };
      return R_dot;
    };
    // need to catch this by reference
    std::vector<Eigen::Triplet<double>> hessian_triplets;
    hessian_triplets.reserve(measurement_size * 4);
    // lambda
    auto get_linear_system_entry = [&](const AtomicMeasurement& z) {
      const auto from_idx = z.from_state_idx;
      const auto to_idx = z.to_state_idx;
      const auto& from = _state[from_idx].rotation;
      const auto& to = _state[to_idx].rotation;
      Eigen::Vector4d h_ij = (from.inverse() * to).toRotationMatrix().reshaped();
      Eigen::Vector4d e = h_ij - z.rotation.toRotationMatrix().reshaped();
      Eigen::Vector4d Ji = (rotation_derived(from.angle()).transpose() * to.toRotationMatrix()).reshaped();
      Eigen::Vector4d Jj = (from.inverse().toRotationMatrix() * rotation_derived(to.angle())).reshaped();
      LinearSystemEntry entry(system_size);
      auto& [g, chi_square] = entry;
      hessian_triplets.emplace_back(from_idx, from_idx, Ji.transpose() * Ji);
      hessian_triplets.emplace_back(to_idx, to_idx, Jj.transpose() * Jj);
      hessian_triplets.emplace_back(from_idx, to_idx, Ji.transpose() * Jj);
      hessian_triplets.emplace_back(to_idx, from_idx, Ji.transpose() * Jj);
      g(from_idx) = -Ji.transpose() * e;
      g(to_idx) = -Jj.transpose() * e;
      chi_square += e.squaredNorm();
      return entry;
    };

    pb_utils::Timer lsq_timer("LSQ Optimization");
    Eigen::VectorXd dx = Eigen::VectorXd::Zero(system_size);
    Eigen::SparseMatrix<double> H(system_size, system_size);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sp_solver;
    /* Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> sp_solver; */
    bool pattern_analyzed = false;
    int iter = 0;
    while (iter < _max_iter) {
      lsq_timer.tick();
      LinearSystemEntry init(system_size);
      // clear the elements, but keep the reserved space
      hessian_triplets.clear();
      /* same for the eigen sparse matrix */
      H.setZero();
      /* do the fishy lambda reduce */
      const auto& [g, chi_square] = std::transform_reduce(_measurement.cbegin(), _measurement.cend(), init,
                                                          std::plus<>(), get_linear_system_entry);
      H.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end());
      const auto& H_block = H.bottomRightCorner(system_size - 1, system_size - 1);
      const auto& g_block = g.tail(system_size - 1);
      if (!pattern_analyzed) {
        /* split the compute step since our sparsity structure remains the same accross iterations */
        sp_solver.analyzePattern(H_block);
        pattern_analyzed = true;
      }
      sp_solver.factorize(H_block);
      dx.tail(system_size - 1) = sp_solver.solve(g_block);
      _state.box_plus(dx);
      auto dx_sqnorm = dx.squaredNorm();

      std::cout << "Iter: " << iter << "/" << _max_iter - 1 << " and chi_squared = " << chi_square
                << " and delta_x (sq. norm) = " << dx_sqnorm << ", took " << lsq_timer.tock() << "s\n";
      if (verbose) {
        std::cout << "g is\n"
                  << g.transpose()
                  /* << "\nHessian is\n"  << H_block << "\nand det is " << H_block.determinant()  */
                  << " \n";
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
