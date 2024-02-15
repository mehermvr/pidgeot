#pragma once

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <numeric>

#include "state.h"
namespace {
struct LinearSystemEntry {
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
  Eigen::Matrix4d H{Eigen::Matrix4d::Zero()};
  Eigen::Vector4d g{Eigen::Vector4d::Zero()};
  double chi_square{0.0};
};

} // namespace
namespace rotsync {
inline auto chi_square(const Eigen::Vector4d& error) { return error.squaredNorm(); };
class Solver {
private:
  State _state;
  Measurement _measurement;
  int _max_iter;
  double _chi_square_thresh;
  int _iter{0};

public:
  Solver(int max_iter, const State& initial_state, const Measurement& measurement, double chi_square_thresh = 1e-8)
      : _state(initial_state), _measurement(measurement), _max_iter(max_iter), _chi_square_thresh(chi_square_thresh) {}

  auto solve(bool verbose = false) {
    auto rotation_derived = [](const double angle) {
      Eigen::Matrix2d R_dot{{-std::sin(angle), -std::cos(angle)}, {std::cos(angle), -std::sin(angle)}};
      return R_dot;
    };
    auto get_linear_system_entry = [&](const int index) {
      const auto& from = _state[index];
      const auto to_index = index < 3 ? index + 1 : 0;
      const auto& to = _state[to_index];
      const auto& z = _measurement(index);
      Eigen::Vector4d h_ij = (from.inverse() * to).toRotationMatrix().reshaped();
      Eigen::Vector4d e = h_ij - Eigen::Rotation2Dd(z).toRotationMatrix().reshaped();
      Eigen::Vector4d Ji = (rotation_derived(from.angle()).transpose() * to.toRotationMatrix()).reshaped();
      Eigen::Vector4d Jj = (from.inverse().toRotationMatrix() * rotation_derived(to.angle())).reshaped();
      LinearSystemEntry entry;
      auto& [H, g, chi_square] = entry;
      H(index, index) = Ji.transpose() * Ji;
      H(to_index, to_index) = Jj.transpose() * Jj;
      H(index, to_index) = Ji.transpose() * Jj;
      H(to_index, index) = H(index, to_index);
      g(index) = -Ji.transpose() * e;
      g(to_index) = -Jj.transpose() * e;
      chi_square += e.squaredNorm();
      return entry;
    };
    auto indices = std::views::iota(0, 4);
    Eigen::Vector4d dx = Eigen::Vector4d::Zero();
    while (_iter < _max_iter) {
      LinearSystemEntry init;
      const auto& [H, g, chi_square] =
          std::transform_reduce(indices.cbegin(), indices.cend(), init, std::plus<>(), get_linear_system_entry);
      dx.tail<3>() = H.bottomRightCorner<3, 3>().ldlt().solve(g.tail<3>());
      _state.box_plus(dx);

      if (verbose) {
        std::cout << "Iter: " << _iter << " and chi_squared = " << chi_square << "\n";
        std::cout << "g is\n"
                  << g.transpose() << "\nHessian is\n"
                  << H.bottomRightCorner<3, 3>() << "\nand det is " << H.bottomRightCorner<3, 3>().determinant()
                  << " \n";
      }
      if (chi_square < _chi_square_thresh) {
        break;
      }
      _iter++;
    }
    return _state;
  }
};
}; // namespace rotsync
