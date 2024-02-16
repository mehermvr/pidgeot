#pragma once

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <numeric>
#include <ranges>

#include "measurement.h"
#include "state.h"

namespace {
struct LinearSystemEntry {
  Eigen::MatrixXd H;
  Eigen::VectorXd g;
  double chi_square{0.0};

  LinearSystemEntry(const size_t size) : H(Eigen::MatrixXd::Zero(size, size)), g(Eigen::VectorXd::Zero(size)) {}

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

public:
  Solver(int max_iter, const State& initial_state, const Measurement& measurement, double chi_square_thresh = 1e-8)
      : _state(initial_state), _measurement(measurement), _max_iter(max_iter), _chi_square_thresh(chi_square_thresh) {}

  auto solve(bool verbose = false) {
    const auto system_size = _measurement.size();
    auto rotation_derived = [](const double angle) {
      Eigen::Matrix2d R_dot{
          {-std::sin(angle), -std::cos(angle)},
          { std::cos(angle), -std::sin(angle)}
      };
      return R_dot;
    };
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
    int iter = 0;
    while (iter < _max_iter) {
      LinearSystemEntry init(system_size);
      const auto& [H, g, chi_square] = std::transform_reduce(_measurement.cbegin(), _measurement.cend(), init,
                                                             std::plus<>(), get_linear_system_entry);
      /* this flops because tail block does not seem assignable */
      /* dx.tail(system_size - 1) = H.bottomRightCorner(system_size -1,
       * system_size-1).ldlt().solve(g.tail(system_size-1)); */
      /* so this is a super ugly hack with a copy, because eigen does not support inserts */
      Eigen::VectorXd dx(system_size);
      dx << 0.0, H.bottomRightCorner(system_size - 1, system_size - 1).ldlt().solve(g.tail(system_size - 1));
      _state.box_plus(dx);

      if (verbose) {
        std::cout << "Iter: " << iter << " and chi_squared = " << chi_square << "\n";
        std::cout << "g is\n"
                  << g.transpose() << "\nHessian is\n"
                  << H.bottomRightCorner(system_size - 1, system_size - 1) << "\nand det is "
                  << H.bottomRightCorner(system_size - 1, system_size - 1).determinant() << " \n";
      }
      if (chi_square < _chi_square_thresh) {
        break;
      }
      iter++;
    }
    return _state;
  }
};
}; // namespace pidgeot
