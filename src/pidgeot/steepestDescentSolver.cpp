#include "steepestDescentSolver.h"
#include "pb_utils/timer.h"

namespace pidgeot {

double SteepestDescentSolver::calculate_step_size(const LinearSystem& linear_system) {
  const auto& g = linear_system.g;
  const auto& H = linear_system.H;
  const double alpha = g.squaredNorm() / (g.transpose() * H * g);
  return alpha;
}

Eigen::VectorXd SteepestDescentSolver::solve(const LinearSystem& linear_system) {
  const double alpha = calculate_step_size(linear_system);
  Eigen::VectorXd sd_step = -alpha * linear_system.g;
  return sd_step;
};
// store the step size in alpha. useful when alpha is needed and to prevent double compute
Eigen::VectorXd SteepestDescentSolver::solve(const LinearSystem& linear_system, double& alpha) {
  alpha = calculate_step_size(linear_system);
  Eigen::VectorXd sd_step = -alpha * linear_system.g;
  return sd_step;
}

State SteepestDescentSolver::solve(bool verbose) {
  pb_utils::Timer lsq_timer("Steepest Descent Optimization");
  const long system_size = pb_utils::saturate_cast<long>(_x.size());
  const long measurement_size = pb_utils::saturate_cast<long>(_measurement.size());
  LinearSystem linear_system(system_size, measurement_size);

  int iter = 0;
  while (iter < _max_iter) {
    if (verbose) {
      lsq_timer.tick();
    }
    linear_system.build(_x, _measurement);
    // use the static overload
    const Eigen::VectorXd dx = solve(linear_system);

    _x.box_plus(dx);
    auto dx_sqnorm = dx.squaredNorm();
    std::cout << "Iter: " << iter << "/" << _max_iter - 1 << " and chi_squared = " << linear_system.chi_square
              << " and delta_x (sq. norm) = " << dx_sqnorm << ", took " << lsq_timer.tock() << "s\n";
    if (linear_system.chi_square < _chi_square_thresh || dx_sqnorm < _dx_sqnorm_thresh) {
      break;
    }
    iter++;
  }
  return _x;
}
} // namespace pidgeot
