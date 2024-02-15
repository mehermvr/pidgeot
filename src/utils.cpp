#include "utils.h"

namespace lutils {
double rad2deg(const double radians) { return radians * 180 / utils::PI; }
double deg2rad(const double degrees) { return degrees * utils::PI / 180; }
} // namespace lutils
