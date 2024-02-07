#include "main.h"

#include <Eigen/src/Geometry/Rotation2D.h>

/* #include <iostream> */
#include <numbers>

auto rotsync::calculate_error_and_jacobian(const State& state, const Measurement& measurement) {
    Jacobian jacobian;
    /* calc error */
    auto wrap_index = [](int idx) { return idx % 4; };
    auto error_scalar = [](const Eigen::Rotation2Df& state_1, const Eigen::Rotation2Df& state_2,
                           const float measurement) {
        auto prediction = state_1.toRotationMatrix().transpose() * state_2.toRotationMatrix();
        auto measurement_so2 = exponential_map(measurement).toRotationMatrix();
        return (Eigen::Matrix3f::Identity() - prediction.transpose() * measurement_so2).trace();
    };
    return std::tuple(error, jacobian);
}

/* int main(int argc, char *argv[]) { */
int main() {
    auto c_pi = std::numbers::pi_v<float>;
    rotsync::Measurement measurement{c_pi, c_pi, c_pi, c_pi};

    rotsync::State initial_state{};

    initial_state.printStateAngles();

    return 0;
}
