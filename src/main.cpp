#include "main.h"

/* #include <iostream> */
#include <numbers>

auto rotsync::calculate_error_and_jacobian(const State& state, const Measurement& measurement) {
    /* calc error */
    auto error_scalar = [](const Eigen::Rotation2Df& state_1, const Eigen::Rotation2Df& state_2,
                           const float relative_rot) {
        auto prediction = state_1.toRotationMatrix().transpose() * state_2.toRotationMatrix();
        auto measurement_so2 = exponential_map(relative_rot).toRotationMatrix();
        return (Eigen::Matrix2f::Identity() - prediction.transpose() * measurement_so2).trace();
    };
    Eigen::Vector4f error;
    Jacobian jacobian = Jacobian::Zero();
    for (int i = 0; i < measurement.size(); i++) {
        error[i] = error_scalar(state[i], state[(i + 1) % 4], measurement[i]);
        jacobian(i, i) = 1;  // needs to be the full element which i dont have here
        jacobian(i, (i + 1) % 4) = -1;
    }
    return std::tuple(error, jacobian);
}

/* int main(int argc, char *argv[]) { */
int main() {
    auto c_pi = std::numbers::pi_v<float>;
    rotsync::Measurement measurement{c_pi, c_pi, c_pi, c_pi};

    rotsync::State initial_state{};

    initial_state.printStateAngles();

    auto [error, jacobian] = rotsync::calculate_error_and_jacobian(initial_state, measurement);

    std::cout << "error is " << error << "\n";
    std::cout << "jacobian is " << jacobian << "\n";
    Eigen::Matrix4f hessian = jacobian.transpose() * jacobian;
    Eigen::Vector4f g = -jacobian.transpose() * error;
    auto delta_x = hessian.llt().solve(g);
    std::cout << delta_x;
    return 0;
}
