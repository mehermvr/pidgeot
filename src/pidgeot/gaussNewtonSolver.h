#pragma once

#include "measurement.h"
#include "state.h"
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <pb_utils/numbers.h>
#include <pb_utils/timer.h>

namespace pidgeot {
class GaussNewtonSolver {
private:
  int _max_iter;
  State _state;
  Measurement _measurement;
  double _chi_square_thresh;
  double _dx_sqnorm_thresh;

public:
  GaussNewtonSolver(int max_iter,
                    const State& initial_state,
                    const Measurement& measurement,
                    double chi_square_thresh = 1e-8,
                    double dx_sqnorm_thresh = 1e-16)
      : _state(initial_state),
        _measurement(measurement),
        _max_iter(max_iter),
        _chi_square_thresh(chi_square_thresh),
        _dx_sqnorm_thresh(dx_sqnorm_thresh) {}
  State solve(bool verbose);
};
}; // namespace pidgeot
