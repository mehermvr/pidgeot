#include "measurement.h"

namespace pidgeot {
// Overload << for Measurement
std::ostream& operator<<(std::ostream& os, const Measurement& measurement) {
  os << "Measurement:";
  for (const auto& atomic : measurement) {
    os << " { " << atomic << " },";
  }
  os << "\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, const AtomicMeasurement& measurement) {
  os << measurement.from_state_idx << " -> " << measurement.to_state_idx << ", Angle: " << measurement.rotation.angle();
  return os;
}
} // namespace pidgeot
