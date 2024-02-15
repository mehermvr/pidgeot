#include <pb_utils/numbers.h>
#include <pb_utils/timer.h>

#include <CLI/CLI.hpp>
#include <iostream>

#include "measurement.h"
/* #include "solver.h" */
#include "state.h"

int main(int argc, char* argv[]) {
  /* CLI start */
  CLI::App app{"Least square kacken"};
  int max_iter = 1000;
  app.add_option("max_iter", max_iter, "Max iterations of Least Squares");
  /* bool use_analytic_jacobian = true; */
  /* app.add_flag("--analytic,!--numeric", use_analytic_jacobian, */
  /* "Use analytic jacobian or numeric."); */
  bool verbose = false;
  app.add_flag("--verbose,-v", verbose, "Additional debug info");
  /* random initial guess that is somewhat close */
  /* std::vector<double> initial_state_vector{0, utils::PI / 1.4, 0.8 * utils::PI, 2 * utils::PI}; */
  /* app.add_option("--initial_guess", initial_state_vector, "Initial state vector. seperate with spaces")->expected(4);
   */
  CLI11_PARSE(app, argc, argv);
  /* CLI done*/

  pigeotto::Measurement measurement{
      {0, 1, utils::PI / 2},
      {1, 2, utils::PI / 2},
      {2, 3, utils::PI / 2},
      {3, 0, utils::PI / 2}
  };
  std::cout << measurement;

  pigeotto::State initial_state({0, 0, 0, 0.1});
  std::cout << initial_state;
  /* std::cout << "Initial state is " << initial_state << "\n"; */

  /* utils::Timer timer("LSQ optimization"); */
  /* pigeotto::Solver solver(max_iter, initial_state, measurement); */
  /* auto final_state = solver.solve(verbose); */

  /* std::cout << "Final state is " << final_state << "\n"; */
  return 0;
}
