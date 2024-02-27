#pragma once
#include "linearSystem.h"
#include "measurement.h"
#include <Eigen/Core>

namespace pidgeot {
class SteepestDescentSolver {
private:
  int _max_iter;
  State _x;
  Measurement _measurement;
  double _chi_square_thresh;
  double _dx_sqnorm_thresh;

public:
  SteepestDescentSolver(int max_iter,
                        const State& initial_state,
                        const Measurement& measurement,
                        double chi_square_thresh = 1e-16,
                        double dx_sqnorm_thresh = 1e-20)
      : _x(initial_state),
        _measurement(measurement),
        _max_iter(max_iter),
        _chi_square_thresh(chi_square_thresh),
        _dx_sqnorm_thresh(dx_sqnorm_thresh) {}

  static double calculate_step_size(const LinearSystem& linear_system);
  static Eigen::VectorXd solve(const LinearSystem& linear_system);
  // step size is stored in alpha
  static Eigen::VectorXd solve(const LinearSystem& linear_system, double& alpha);
  State solve(bool verbose = false);
};
} // namespace pidgeot
