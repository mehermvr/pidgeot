#pragma once
#include <Eigen/Geometry>
#include <pb_utils/numbers.h>
#include <vector>

namespace pidgeot {

struct AtomicMeasurement {
  int from_state_idx;
  int to_state_idx;
  Eigen::Rotation2Dd rotation;

  AtomicMeasurement(const int from, const int to, const double angle)
      : from_state_idx(from), to_state_idx(to), rotation(Eigen::Rotation2Dd(angle)) {}

  friend std::ostream& operator<<(std::ostream& os, const AtomicMeasurement& measurement) {
    os << measurement.from_state_idx << " -> " << measurement.to_state_idx
       << ", Angle: " << measurement.rotation.angle() << "\u00B0";
    return os;
  }
};

using Measurement = std::vector<AtomicMeasurement>;
// Overload << for Measurement
std::ostream& operator<<(std::ostream& os, const Measurement& measurement);
} // namespace pidgeot
