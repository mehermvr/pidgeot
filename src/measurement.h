#pragma once
#include "utils.h"
#include <Eigen/Geometry>
#include <vector>

namespace pigeotto {

struct AtomicMeasurement {
  int from_state_idx;
  int to_state_idx;
  Eigen::Rotation2Dd transform;

  AtomicMeasurement(const int from, const int to, const double angle)
      : from_state_idx(from), to_state_idx(to), transform(Eigen::Rotation2Dd(angle)) {}

  friend std::ostream& operator<<(std::ostream& os, const AtomicMeasurement& measurement) {
    os << measurement.from_state_idx << " -> " << measurement.to_state_idx
       << ", Angle: " << lutils::rad2deg(measurement.transform.angle()) << "\u00B0";
    return os;
  }
};

using Measurement = std::vector<AtomicMeasurement>;
// Overload << for Measurement
std::ostream& operator<<(std::ostream& os, const Measurement& measurement);
} // namespace pigeotto
