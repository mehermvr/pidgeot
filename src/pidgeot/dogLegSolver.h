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
/* see Methods for Non-Linear Least Squares Problems, Madsen et al. 2004 */
/* Termination Criteria - Powell's Dog Leg Solver
 * 1. inifinity_norm(g) <= eps_1 => gradient is close to 0, we are at a local minimum
 * 2. (dx).norm() <= eps_2 => absolute change in x is small. NOTE this is different from madsen's termination criterion
 * which is based on relative change
 * 3. iter >= iter_max => prevent inifinite iterations
 * 4. inifinity_norm(error) <= eps_3 => error is close to 0. We do a threshold on chi_square instead since this is
 * easier in the current implementation
 * 5. trust_radius <= eps_2(x.norm + eps_2) => because if this is true, condition 2 is gauranteed. NOTE this is not
 * implemented since this also involves x.norm like criterion 2 which is a pain. its gauranteed, let the next iteration
 * handle it
 */
class DogLegSolver {
public:
  using SparseSolver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;
  /*using SparseSolver = Eigen::ConjugateGradient<Eigen::SparseMatrix<double>>; */

  /*
   * param eps_1: if max element in gradient at an iteration is less than eps_1, terminate
   * param eps_2: if obsolute change in x across an interation (dx) is less than eps_2, terminate. this is different
   * from madsen, this the absolute is easier to implement than relative. (state_norm is a pain)
   * param eps_3: if chi_square in an iteration is less than eps_3, terminate. this is different from Madsen
   */
  DogLegSolver(int max_iter,
               const State& initial_state,
               const Measurement& measurement,
               double trust_radius = 1e4,
               double eps_1 = 1e-4,  // gradient infinity norm threshold
               double eps_2 = 1e-12, // dx norm absolute threshold
               double eps_3 = 1e-8   // chi_square threshold
               )
      : _x(initial_state),
        _measurement(measurement),
        _max_iter(max_iter),
        _trust_radius(trust_radius),
        _eps_1(eps_1),
        _eps_2(eps_2),
        _eps_3(eps_3) {}

  /* returns the final state */
  State solve(bool verbose);

private:
  void _calculate_trust_diagonal(const LinearSystem& linear_system);

  Eigen::VectorXd _trust_scale_vector(const Eigen::VectorXd& vec) const;
  Eigen::VectorXd _trust_unscale_vector(const Eigen::VectorXd& scaled_vec) const;

  Eigen::VectorXd _calculate_gauss_newton_step(const LinearSystem& linear_system);

  std::pair<double, Eigen::VectorXd> _calculate_cauchy_step(const LinearSystem& linear_system);

  /* returns dl_step, the dl switch case, and the beta value. beta is only meaningful in case 3 */
  /* we keep track of which dl case was the iteration, so we can calculate the linear model diff L(0) - L(dl_dx)
   * without keeping track of the Jacobian */
  std::tuple<Eigen::VectorXd, int, double> _calculate_dogleg_step(const auto& sd_step, const auto& gn_step);

  /* gain ratio between the non linear model and the linear model approximation */
  double _calculate_gain_ratio(const State& x_new,
                               const LinearSystem& linear_system,
                               const int dl_case,
                               const double alpha,
                               const double beta) const;

  int _max_iter;
  State _x;
  Measurement _measurement;
  double _eps_1; // gradient infinity norm threshold
  double _eps_2; // dx norm absolute threshold
  double _eps_3; // chi_square threshold
  double _trust_radius;
  Eigen::VectorXd _trust_diagonal;
  SparseSolver _sparse_solver;
  bool _sparse_pattern_analyzed{false};
};
}; // namespace pidgeot
