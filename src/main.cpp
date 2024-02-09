#include "main.h"

#include <CLI/App.hpp>
#include <CLI/CLI.hpp>
/* #include <iostream> */
#include <numbers>

Eigen::Vector4d calculate_error(const rotsync::State& state,
                                const rotsync::Measurement& measurement) {
    auto calc_error_scalar = [](const Eigen::Rotation2Dd& state_1,
                                const Eigen::Rotation2Dd& state_2, const double relative_rot) {
        Eigen::Matrix2d prediction =
                state_1.toRotationMatrix().transpose() * state_2.toRotationMatrix();
        Eigen::Matrix2d measurement_so2 = Eigen::Rotation2Dd(relative_rot).toRotationMatrix();
        return (Eigen::Matrix2d::Identity() - prediction.transpose() * measurement_so2).trace();
    };
    Eigen::Vector4d error = Eigen::Vector4d::Zero();
    for (int i = 0; i < error.size(); i++) {
        error[i] = calc_error_scalar(state[i], state[(i + 1) % 4], measurement[i]);
    }
    return error;
}
Eigen::Matrix4d calculate_numerical_jacobian(const rotsync::State& state,
                                             const rotsync::Measurement& measurement) {
    Eigen::Matrix4d jacobian = Eigen::Matrix4d::Zero();
    auto initial_error = calculate_error(state, measurement);
    double epsilon = 1e-8;
    for (int idx = 0; idx < jacobian.cols(); idx++) {
        Eigen::Vector4d perturbation = Eigen::Vector4d::Zero();
        perturbation[idx] += epsilon;
        auto plussed_state = state.box_plus(perturbation);
        auto perturbed_error = calculate_error(plussed_state, measurement);
        auto delta_error = perturbed_error - initial_error;
        jacobian.col(idx) = delta_error / epsilon;
    }
    return jacobian;
}
Eigen::Matrix4d calculate_analytical_jacobian(const rotsync::State& state,
                                              const rotsync::Measurement& measurement) {
    Eigen::Matrix4d jacobian = Eigen::Matrix4d::Zero();
    auto wrap_index = [](int idx) { return idx % 4; };
    for (int idx = 0; idx < jacobian.cols(); idx++) {
        jacobian(idx, idx) = 2 * std::sin(state[idx].angle() - state[wrap_index(idx + 1)].angle() +
                                          measurement[idx]);
        jacobian(idx, wrap_index(idx + 1)) =
                -2 * std::sin(state[idx].angle() - state[wrap_index(idx + 1)].angle() +
                              measurement[idx]);
    }
    return jacobian;
}
auto rotsync::calculate_error_and_jacobian(const State& state,
                                           const Measurement& measurement,
                                           bool use_analytic_jacobian) {
    /* calc error */
    Eigen::Vector4d error = calculate_error(state, measurement);
    Jacobian jacobian;
    if (use_analytic_jacobian) {
        jacobian = calculate_analytical_jacobian(state, measurement);
    } else {
        jacobian = calculate_numerical_jacobian(state, measurement);
    }
    return std::tuple(error, jacobian);
}

int main(int argc, char* argv[]) {
    /* CLI start */
    CLI::App app{"Least square kacken"};
    int max_iter = 1000;
    app.add_option("max_iter", max_iter, "Max iterations of Least Squares");
    bool use_analytic_jacobian = true;
    app.add_flag("--analytic,!--numeric", use_analytic_jacobian,
                 "Use analytic jacobian or numeric.");

    CLI11_PARSE(app, argc, argv);
    /* CLI done*/

    auto c_pi = std::numbers::pi_v<double>;
    rotsync::Measurement measurement{c_pi / 2, c_pi / 2, c_pi / 2, c_pi / 2};

    rotsync::State initial_state{};
    /* rotsync::State initial_state{{0, c_pi / 2, c_pi, 3 * c_pi / 2}}; */
    std::cout << "Initial state is " << initial_state << "\n";

    auto state = initial_state;
    int iter = 0;
    auto chi_square = [](const Eigen::Vector4d& error) { return error.squaredNorm(); };
    double chi_square_thresh = 1e-3;
    while (iter < max_iter) {
        auto [error, jacobian] =
                rotsync::calculate_error_and_jacobian(state, measurement, use_analytic_jacobian);

        std::cout << "error is " << error << "\n";
        std::cout << "jacobian is\n"
                  << jacobian << "\nand determinant is " << jacobian.determinant() << " \n";
        Eigen::Matrix<double, 4, 3> jacobian_with_fixed_dof = jacobian.block<4, 3>(0, 1);
        std::cout << "jacobian with fixed dof is\n" << jacobian_with_fixed_dof << "\n";
        Eigen::Matrix3d hessian = jacobian_with_fixed_dof.transpose() * jacobian_with_fixed_dof;
        Eigen::Vector3d rhs = -jacobian_with_fixed_dof.transpose() * error;
        Eigen::Vector4d delta_x = Eigen::Vector4d::Zero();
        /* first measurement is fixed */
        delta_x.tail(3) = hessian.llt().solve(rhs);
        std::cout << "hessian \n"
                  << hessian << "\nand determinant is " << hessian.determinant() << " \n";
        std::cout << "rhs " << rhs << "\n";
        /* std::cout << "delta_x is " << delta_x << "\n"; */
        /* std::cout << "error is " << error << "\n"; */
        state = state.box_plus(delta_x);
        if (chi_square(error) < chi_square_thresh) {
            /* this placed here means one extra iteration though */
            std::cout << "error is " << error << "\n";
            std::cout << "jacobian is\n"
                      << jacobian << "\nand determinant is " << jacobian.determinant() << " \n";
            std::cout << "hessian is\n"
                      << hessian << "\nand det is " << hessian.determinant() << " \n";
            break;
        }
        iter++;
    }
    std::cout << "Final state is " << state << "\n";
    return 0;
}
