#pragma once

#include "linearSystem.h"
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
// use with a linear_system whose sparse structure doesn't change. otherwise needs to be adapted
class GaussNewtonSolver {
public:
  using SparseSolver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;
  /*using SparseSolver = Eigen::ConjugateGradient<Eigen::SparseMatrix<double>>; */

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

  static Eigen::VectorXd
  solve(const LinearSystem& linear_system, SparseSolver& sparse_solver, bool& sparse_pattern_analyzed);

private:
  int _max_iter;
  State _state;
  Measurement _measurement;
  double _chi_square_thresh;
  double _dx_sqnorm_thresh;
  bool _sp_pattern_analyzed{false};
  SparseSolver _sp_solver;
};
;
}; // namespace pidgeot
