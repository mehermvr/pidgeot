#pragma once
#include <pb_utils/numbers.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <iostream>
#include <ranges>

namespace rotsync {
using Perturbation = Eigen::Vector4d;
using Measurement = Eigen::Vector4d;
using Jacobian = Eigen::Matrix4d;

inline Eigen::Rotation2Dd exponential_map(const double angle) { return Eigen::Rotation2Dd(angle); }

class State {
public:
  explicit State(size_t size) : _states_so2(size, Eigen::Rotation2Dd(0.0)) {}

  /* Construct with a container of angles. assumes radians input. */
  template <typename T>
    requires std::ranges::input_range<T>
  explicit State(const T& angles) {
    std::ranges::for_each(angles, [this](const auto angle) { _states_so2.emplace_back(angle); });
  }

  /* convert the so2 state vector to an eiger euler angle vector and return */
  Eigen::Vector4d to_angles() const {
    Eigen::Vector4d angles;
    std::ranges::transform(_states_so2, angles.begin(),
                           [](const auto& rot_mat) { return rot_mat.angle(); });
    return angles;
  }

  const auto& as_vector() const { return _states_so2; }

  /* box plus with the perturbation R4 vector, return a new state*/
  void box_plus(const Perturbation& perturbation) {
    std::ranges::transform(_states_so2, perturbation, _states_so2.begin(),
                           [](const Eigen::Rotation2Dd& state, const double pert_angle) {
                             return exponential_map(pert_angle) * state;
                           });
  }
  /* const reference to an element of the so2 state vector */
  const Eigen::Rotation2Dd& operator[](const std::size_t idx) const {
    assert(idx < 4 && "Index must be less than 4");
    return _states_so2.at(idx);
  }
  // Overload << operator for printing
  friend std::ostream& operator<<(std::ostream& ostream, const State& state) {
    ostream << "State {";
    const auto& angles = state.to_angles();
    for (const auto& angle : state.to_angles()) {
      ostream << angle * 180 / utils::PI << "\u00B0,";
    }
    ostream << "}";
    return ostream;
  }

private:
  std::vector<Eigen::Rotation2Dd> _states_so2;
};

} // namespace rotsync
