#include <algorithm>
#include <array>
#include <pb_utils/numbers.h>
#include <pb_utils/timer.h>

#include "pidgeot/dogLegSolver.h"
#include "pidgeot/gaussNewtonSolver.h"
#include "pidgeot/measurement.h"
#include "pidgeot/state.h"
#include <CLI/CLI.hpp>
#include <iostream>
#include <numbers>
#include <random>

int main(int argc, char* argv[]) {
  pb_utils::Timer timer("Program execution");
  /* CLI start */
  CLI::App app{"Pidgeot"};
  int seed = 42;
  app.add_option("--seed,-s", seed, "set the seed");
  int max_iter = 1000;
  app.add_option("max_iter", max_iter, "Max iterations of Least Squares");
  int state_length = 100;
  app.add_option("--state-length,-l", state_length, "Length of the state vector");
  int max_loop_closures = 10;
  app.add_option("--max_lc,-c", max_loop_closures, "Maximum number of loop closures");
  /* unused. but keep for later */
  bool use_analytic_jacobian = true;
  app.add_flag("--analytic,!--numeric", use_analytic_jacobian, "Use analytic jacobian or numeric.");
  bool verbose = false;
  app.add_flag("--verbose,-v", verbose, "Additional debug info");
  CLI11_PARSE(app, argc, argv);
  /* CLI done*/

  using std::numbers::pi;

  /* can get an rng with implementation defined seed state mt_engine{std::random_device} or can use mt_engine.seed(seed)
   * for determinism. but use a seed seq to get good entropy, and still be deterministic */
  std::seed_seq seed_seq{seed};
  std::mt19937 mt_engine{seed_seq};

  /* generate a length of random states */
  std::uniform_real_distribution<double> angle_range{0, 2 * pi};
  std::vector<double> gt_state_angles(state_length);
  /* first state is 0, since we fix it later */
  std::ranges::generate(std::next(gt_state_angles.begin()), gt_state_angles.end(),
                        [&]() { return angle_range(mt_engine); });
  std::cout << "Ground truth state length is " << gt_state_angles.size() << "\n";

  /* generate measurements */
  pidgeot::Measurement measurement;
  /* precalculate the expected length to reserve the space before emplace backing */
  const auto actual_max_pairs = (state_length * (state_length - 1) / 2) - (state_length - 1);
  if (actual_max_pairs < max_loop_closures) {
    std::cout << "Max_loop_closures " << max_loop_closures << " > " << actual_max_pairs
              << " actual_max_pairs possible\nSetting it to limit\n";
    max_loop_closures = actual_max_pairs;
  }
  const auto loop_closures = std::uniform_int_distribution<int>{max_loop_closures / 3, max_loop_closures}(mt_engine);
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
  std::cout << "Measurment vector length is " << measurement.size() << "\n";

  pidgeot::State initial_state(state_length);
  std::cout << "Initial state length is " << initial_state.size() << "\n";

  pidgeot::DogLegSolver solver(max_iter, initial_state, measurement);
  auto final_state = solver.solve(verbose);

  double error = 0;
  std::ranges::for_each(std::views::enumerate(gt_state_angles), [&](const auto& enum_elem) {
    const auto& [idx, gt_angle] = enum_elem;
    const auto pred_angle = final_state[idx].rotation.smallestPositiveAngle();
    const auto diff = pred_angle - gt_angle;
    error += abs(diff);
  });
  std::cout << "mean abs radian error " << error / state_length << " = " << pb_utils::rad2deg(error / state_length)
            << " deg\n";
  return 0;
}
