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

} // namespace pidgeot
