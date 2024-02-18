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
} // namespace pidgeot
