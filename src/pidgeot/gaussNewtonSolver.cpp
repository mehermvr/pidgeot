#include "gaussNewtonSolver.h"

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
    dx = solve(linear_system, _sp_solver, _sp_pattern_analyzed);
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
/* calculate the gauss newton step. if the sparse pattern was not analyzed, sparse_pattern_analyzed is set to true */
Eigen::VectorXd GaussNewtonSolver::solve(const LinearSystem& linear_system,
                                         SparseSolver& sparse_solver,
                                         bool& sparse_pattern_analyzed) {
  const auto& H = linear_system.H;
  const auto& g = linear_system.g;
  /* the system is (H:= J.T*J) * dx = -g */
  const auto system_size = H.outerSize();
  if (!sparse_pattern_analyzed) {
    /* split the compute step since our sparsity structure remains the same accross iterations */
    sparse_solver.analyzePattern(H);
    sparse_pattern_analyzed = true;
  }
  sparse_solver.factorize(H);
  Eigen::VectorXd dx = Eigen::VectorXd::Zero(system_size);
  dx = sparse_solver.solve(-1 * g);
  return dx;
}

} // namespace pidgeot
