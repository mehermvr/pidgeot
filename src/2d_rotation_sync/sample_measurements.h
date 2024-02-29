#pragma once
#include "pidgeot/measurement.h"
#include <vector>

/* splitting stuff to improve compilation speeds */
pidgeot::Measurement sample_measurements(const int seed,
                                         const long state_length,
                                         long max_loop_closures,
                                         const std::vector<double>& gt_state_angles);
