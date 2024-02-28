#include "linearSystem.h"

namespace pidgeot {
/* call before solve */
void LinearSystem::build(const State& x, const Measurement& measurement) {
  clear();
  auto rotation_derived = [](const double angle) {
    Eigen::Matrix2d R_dot{
        {-std::sin(angle), -std::cos(angle)},
        { std::cos(angle), -std::sin(angle)}
    };
    return R_dot;
  };
  auto update_linear_system = [&](const AtomicMeasurement& z) {
    const auto from_idx = z.from_state_idx;
    const auto to_idx = z.to_state_idx;
    const auto& from = x[from_idx];
    const auto& to = x[to_idx];
    /* const auto& [from, to] = z.get_associated_states(x); */
    const Eigen::Vector4d e = z.calculate_error(from, to);
    Eigen::Vector4d Ji =
        (rotation_derived(from.rotation.angle()).transpose() * to.rotation.toRotationMatrix()).reshaped();
    Eigen::Vector4d Jj =
        (from.rotation.inverse().toRotationMatrix() * rotation_derived(to.rotation.angle())).reshaped();
    hessian_triplets.emplace_back(from_idx, from_idx, Ji.transpose() * Ji);
    hessian_triplets.emplace_back(to_idx, to_idx, Jj.transpose() * Jj);
    hessian_triplets.emplace_back(from_idx, to_idx, Ji.transpose() * Jj);
    hessian_triplets.emplace_back(to_idx, from_idx, Ji.transpose() * Jj);
    /* g is the gradient, not the rhs in the normal equations */
    g(from_idx) += Ji.transpose() * e;
    g(to_idx) += Jj.transpose() * e;
    chi_square += e.squaredNorm();
  };
  std::ranges::for_each(measurement, update_linear_system);
  /* std::cout << g << "\n"; */
  convert_triplets();
}
void LinearSystem::clear() {
  hessian_triplets.clear();
  H.setZero();
  g.setZero();
  chi_square = 0.0;
}
void LinearSystem::convert_triplets() { H.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end()); }

/* calculate the gauss newton step */
Eigen::VectorXd LinearSystem::solve() {
  /* the system is (H:= J.T*J) * dx = -g */
  const auto system_size = H.outerSize();
  /* TODO: this is because we fix the first state. there has to be a better way */
  const auto& H_block = H.bottomRightCorner(system_size - 1, system_size - 1);
  const auto& g_block = g.tail(system_size - 1);
  if (!sp_pattern_analyzed) {
    /* split the compute step since our sparsity structure remains the same accross iterations */
    sp_solver.analyzePattern(H_block);
    sp_pattern_analyzed = true;
  }
  sp_solver.factorize(H_block);
  Eigen::VectorXd dx = Eigen::VectorXd::Zero(system_size);
  dx.tail(system_size - 1) = sp_solver.solve(-1 * g_block);
  return dx;
}

} // namespace pidgeot
