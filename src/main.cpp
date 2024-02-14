#include <pb_utils/numbers.h>
#include <pb_utils/timer.h>

#include <CLI/CLI.hpp>
#include <iostream>

#include "solver.h"
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
    /* error at 0 is 0 */
    /* std::array<double, 4> initial_state_array{0, utils::C_PI / 2, utils::C_PI, 3 * utils::C_PI /
     * 2}; */
    /* random initial guess that is somewhat close */
    std::vector<double> initial_state_vector{0, utils::C_PI / 1.4, 0.8 * utils::C_PI,
                                             2 * utils::C_PI};
    app.add_option("--initial_guess", initial_state_vector,
                   "Initial state vector. seperate with spaces")
            ->expected(4);
    CLI11_PARSE(app, argc, argv);
    /* CLI done*/

    rotsync::Measurement measurement{utils::C_PI / 2, utils::C_PI / 2, utils::C_PI / 2,
                                     utils::C_PI / 2};

    rotsync::State initial_state{initial_state_vector};
    std::cout << "Initial state is " << initial_state << "\n";

    utils::Timer timer("LSQ optimization");
    rotsync::Solver solver(initial_state, measurement, max_iter);
    auto final_state = solver.solve(verbose);
    std::cout << "Final state is " << final_state << "\n";
    return 0;
}
