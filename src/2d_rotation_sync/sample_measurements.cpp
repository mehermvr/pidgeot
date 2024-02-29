#include "sample_measurements.h"
#include "pidgeot/measurement.h"
#include <random>

pidgeot::Measurement sample_measurements(const int seed,
                                         const long state_length,
                                         long max_loop_closures,
                                         const std::vector<double>& gt_state_angles) {
  /* generate measurements */
  pidgeot::Measurement measurement;
  /* precalculate the expected length to reserve the space before emplace backing */
  const auto actual_max_pairs = (state_length * (state_length - 1) / 2) - (state_length - 1);
  if (actual_max_pairs < max_loop_closures) {
    std::cout << "Max_loop_closures " << max_loop_closures << " > " << actual_max_pairs
              << " actual_max_pairs possible\nSetting it to limit\n";
    max_loop_closures = actual_max_pairs;
  }

  std::seed_seq seed_seq{seed};
  std::mt19937 mt_engine{seed_seq};
  const auto loop_closures = std::uniform_int_distribution<long>{max_loop_closures / 3, max_loop_closures}(mt_engine);
  std::cout << "Max loop closures: " << max_loop_closures << "\n";
  std::cout << "Sampled loop closures: " << loop_closures << "\n";
  /* state_length - 1 sequential measurements and loop_closures loop closing measurements */
  measurement.reserve(state_length - 1 + loop_closures);

  // reused lamda
  auto get_relative_angle = [&](const int from_idx, const int to_idx) {
    const auto from_angle = gt_state_angles.at(from_idx);
    const double to_angle = gt_state_angles.at(to_idx);
    const double relative_angle =
        (Eigen::Rotation2Dd(from_angle).inverse() * Eigen::Rotation2Dd(to_angle)).smallestPositiveAngle();
    return relative_angle;
  };
  /* generate the sequential measurements */
  std::ranges::for_each(std::views::iota(0, state_length - 1), [&](const int from_idx) {
    const auto to_idx = from_idx + 1;
    const auto relative_angle = get_relative_angle(from_idx, to_idx);
    measurement.emplace_back(from_idx, to_idx, relative_angle);
  });
  /* generate loop closures */
  std::array<int, 2> sampled_indices{};
  const auto idx_range = std::views::iota(0, state_length);
  for (int i = 0; i < loop_closures; i++) {
    std::ranges::sample(idx_range, sampled_indices.begin(), 2, mt_engine);
    const auto [from_idx, to_idx] = sampled_indices;
    const auto relative_angle = get_relative_angle(from_idx, to_idx);
    measurement.emplace_back(from_idx, to_idx, relative_angle);
  }
  return measurement;
}
