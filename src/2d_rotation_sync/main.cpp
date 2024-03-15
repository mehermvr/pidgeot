#include <algorithm>
#include <pb_utils/numbers.h>
#include <pb_utils/timer.h>

#include "pidgeot/dogLegSolver.h"
#include "pidgeot/gaussNewtonSolver.h"
#include "pidgeot/measurement.h"
#include "pidgeot/state.h"
#include "pidgeot/steepestDescentSolver.h"
#include "sample_measurements.h"
#include "sample_states.h"
#include <CLI/CLI.hpp>
#include <iostream>

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

  std::vector<double> gt_state_angles = sample_states(seed, state_length);
  std::cout << "Ground truth state length is " << gt_state_angles.size() << "\n";

  pidgeot::Measurement measurement = sample_measurements(seed, state_length, max_loop_closures, gt_state_angles);
  std::cout << "Measurment vector length is " << measurement.size() << "\n";

  pidgeot::State initial_state(state_length);
  std::cout << "Initial state length is " << initial_state.size() << "\n";
  std::cout << "Fixing state 0\n";
  initial_state.fix_state(0);

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
  /* std::cout << "Final state is " << final_state; */
  return 0;
}
