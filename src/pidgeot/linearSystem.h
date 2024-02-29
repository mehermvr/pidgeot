#pragma once
#include "measurement.h"
#include "state.h"
#include <Eigen/Core>
#include <Eigen/SparseCholesky>

namespace pidgeot {
class LinearSystem {
public:
  /* hessian, but more like J.T * J, ignores the second order terms. if the error at optimum is 0, then H = J.T * J and
   * you get quadratic convergence with GN. see madsen */
  Eigen::SparseMatrix<double> H;
  /* the gradient. the system is (J.T*J)h = -g */
  Eigen::VectorXd g;
  double chi_square{0.0};

private:
  std::vector<Eigen::Triplet<double>> hessian_triplets;

public:
  LinearSystem(const long system_size, const long measurement_size)
      : g(Eigen::VectorXd::Zero(system_size)), H(system_size, system_size) {
    /* allocates memory for hessian triples, the sparse hessian matrices */
    hessian_triplets.reserve(measurement_size * 4);
  }
  /* call before solve */
  void build(const State& x, const Measurement& measurement);
  /* calculate the gauss newton step */
  Eigen::VectorXd solve();

private:
  /* clears all the entries to allow reusing the allocated memory */
  void clear();
  void convert_triplets();
};

} // namespace pidgeot
