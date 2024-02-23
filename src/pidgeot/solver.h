#pragma once

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <numeric>

#include "measurement.h"
#include "state.h"
#include <pb_utils/numbers.h>
#include <pb_utils/timer.h>

namespace {
struct LinearSystemEntry {
  Eigen::MatrixXd H;
  Eigen::VectorXd g;
  double chi_square{0.0};

  explicit LinearSystemEntry(const long size) : H(Eigen::MatrixXd::Zero(size, size)), g(Eigen::VectorXd::Zero(size)) {}

  LinearSystemEntry& operator+=(const LinearSystemEntry& other) {
    this->H += other.H;
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
         double dx_sqnorm_thresh = 1e-8)
      : _state(initial_state),
        _measurement(measurement),
        _max_iter(max_iter),
        _chi_square_thresh(chi_square_thresh),
        _dx_sqnorm_thresh(dx_sqnorm_thresh) {}

  auto solve(bool verbose = false) {
    const long system_size = pb_utils::saturate_cast<long>(_state.size());

    // lambda
    auto rotation_derived = [](const double angle) {
      Eigen::Matrix2d R_dot{
          {-std::sin(angle), -std::cos(angle)},
          { std::cos(angle), -std::sin(angle)}
      };
      return R_dot;
    };
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
      auto& [H, g, chi_square] = entry;
      H(from_idx, from_idx) = Ji.transpose() * Ji;
      H(to_idx, to_idx) = Jj.transpose() * Jj;
      H(from_idx, to_idx) = Ji.transpose() * Jj;
      H(to_idx, from_idx) = H(from_idx, to_idx);
      g(from_idx) = -Ji.transpose() * e;
      g(to_idx) = -Jj.transpose() * e;
      chi_square += e.squaredNorm();
      return entry;
    };

    pb_utils::Timer lsq_timer("LSQ Optimization");
    Eigen::VectorXd dx = Eigen::VectorXd::Zero(system_size);
    int iter = 0;
    while (iter < _max_iter) {
      lsq_timer.tick();
      LinearSystemEntry init(system_size);
      const auto& [H, g, chi_square] = std::transform_reduce(_measurement.cbegin(), _measurement.cend(), init,
                                                             std::plus<>(), get_linear_system_entry);
      const auto& H_block = H.bottomRightCorner(system_size - 1, system_size - 1);
      const auto& g_block = g.tail(system_size - 1);
      dx.tail(system_size - 1) = H_block.ldlt().solve(g_block);
      _state.box_plus(dx);
      auto dx_sqnorm = dx.squaredNorm();

      std::cout << "Iter: " << iter << "/" << _max_iter - 1 << " and chi_squared = " << chi_square
                << " and delta_x (sq. norm) = " << dx_sqnorm << ", took " << lsq_timer.tock() << "s\n";
      if (verbose) {
        std::cout << "g is\n"
                  << g.transpose() << "\nHessian is\n"
                  << H_block << "\nand det is " << H_block.determinant() << " \n";
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
