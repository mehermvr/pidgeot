#pragma once
#include "state.h"
#include <Eigen/Geometry>
#include <pb_utils/numbers.h>
#include <tuple>
#include <vector>
namespace pidgeot {

struct AtomicMeasurement {
  using Rotation = Eigen::Rotation2Dd;
  int from_state_idx;
  int to_state_idx;
  Rotation rotation;

  AtomicMeasurement(const int from, const int to, const double angle)
      : from_state_idx(from), to_state_idx(to), rotation(Eigen::Rotation2Dd(angle)) {}

  friend std::ostream& operator<<(std::ostream& os, const AtomicMeasurement& measurement);

  std::tuple<const AtomicState&, const AtomicState&> get_associated_states(const State& x) const {
    const auto& from = x[from_state_idx];
    const auto& to = x[to_state_idx];
    return std::make_tuple(from, to);
  }
  Eigen::Vector4d calculate_error(const AtomicState& from, const AtomicState& to) const {
    const Eigen::Vector4d h_ij = (from.rotation.inverse() * to.rotation).toRotationMatrix().reshaped();
    Eigen::Vector4d e = h_ij - rotation.toRotationMatrix().reshaped();
    return e;
  };
};

using Measurement = std::vector<AtomicMeasurement>;
// Overload << for Measurement
std::ostream& operator<<(std::ostream& os, const Measurement& measurement);
} // namespace pidgeot
