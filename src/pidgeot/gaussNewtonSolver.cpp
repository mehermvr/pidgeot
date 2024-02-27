#include "gaussNewtonSolver.h"
#include "pidgeot/linearSystem.h"
#include <Eigen/src/Core/Matrix.h>
#include <numeric>

namespace {
struct LinearSystemEntry {
  /* Eigen::MatrixXd H; */
  Eigen::VectorXd g;
  double chi_square{0.0};

  explicit LinearSystemEntry(const long size) : g(Eigen::VectorXd::Zero(size)) {}

  LinearSystemEntry& operator+=(const LinearSystemEntry& other) {
    this->g += other.g;
    this->chi_square += other.chi_square;
    return *this;
  }

  friend LinearSystemEntry operator+(LinearSystemEntry lhs, const LinearSystemEntry& rhs) {
    lhs += rhs;
    return lhs;
  }
};
} // namespace

namespace pidgeot {

State GaussNewtonSolver::solve(bool verbose) {
  pb_utils::Timer lsq_timer("Gauss Newton Optimization");
  if (verbose) {
    lsq_timer.tick();
  }
  const long system_size = pb_utils::saturate_cast<long>(_state.size());
  const long measurement_size = pb_utils::saturate_cast<long>(_measurement.size());

  LinearSystem linear_system(system_size, measurement_size);
  Eigen::VectorXd dx = Eigen::VectorXd::Zero(system_size);
  int iter = 0;
  while (iter < _max_iter) {
    lsq_timer.tick();
    linear_system.build(_state, _measurement);
    dx = linear_system.solve();
    _state.box_plus(dx);
    auto dx_sqnorm = dx.squaredNorm();

    std::cout << "Iter: " << iter << "/" << _max_iter - 1 << " and chi_squared = " << linear_system.chi_square
              << " and delta_x (sq. norm) = " << dx_sqnorm;
    if (verbose) {
      std::cout << ", took " << lsq_timer.tock() << "s";
    }
    std::cout << "\n";
    if (linear_system.chi_square < _chi_square_thresh || dx_sqnorm < _dx_sqnorm_thresh) {
      break;
    }
    iter++;
  }
  return _state;
}

} // namespace pidgeot
