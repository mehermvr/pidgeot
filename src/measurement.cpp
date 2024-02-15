#include "measurement.h"
#include "utils.h"

namespace pigeotto {
// Overload << for Measurement
std::ostream& operator<<(std::ostream& os, const Measurement& measurement) {
  os << "Measurement:";
  for (const auto& atomic : measurement) {
    os << " { " << atomic << " },";
  }
  os << "\n";
  return os;
}
} // namespace pigeotto
