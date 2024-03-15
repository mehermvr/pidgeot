#include "sample_states.h"
#include <algorithm>
#include <numbers>
#include <random>

std::vector<double> sample_states(const int seed, const long state_length) {
  using std::numbers::pi;

  /* can get an rng with implementation defined seed state mt_engine{std::random_device} or can use mt_engine.seed(seed)
   * for determinism. but use a seed seq to get good entropy, and still be deterministic */
  std::seed_seq seed_seq{seed};
  std::mt19937 mt_engine{seed_seq};

  /* generate a length of random states */
  std::uniform_real_distribution<double> angle_range{0, 2 * pi};
  std::vector<double> gt_state_angles(state_length);
  /* first state is 0, since we fix it later and to allow easy eval with initial state
   * of 0 */
  std::ranges::generate(std::next(gt_state_angles.begin()), gt_state_angles.end(),
                        [&]() { return angle_range(mt_engine); });
  return gt_state_angles;
}
