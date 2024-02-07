#pragma once

/* the document */
/*
 * H * delta_x = -g
 * H = J' J
 * g = J' error
 * x' = x + delta_x
 * iterate till convergence
 *
 * state SO2 x4
 * perturbation R4
 * box_plus(state, perturbation) = Exp(perturbation) * state
 * prediction(state, indexe i) = R_(i-1)' * R_i
 * measurement(i)
 * error(state, i) = tr(1 - prediction(state, indexed i)'measurement(i))
 * Jacobian(state, perturbation, measurement, i): R1x4 -> numerical or analytic
 * go to solve lin system
 */

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <iostream>
#include <tuple>

namespace rotsync {
using Perturbation = Eigen::Vector4f;
using Measurement = Eigen::Vector4f;
using Jacobian = Eigen::Matrix4f;

inline Eigen::Rotation2Df exponential_map(const float angle) { return Eigen::Rotation2Df(angle); }

class State {
public:
    /* default constructor, all absolute rotations are 0 rad */
    State() { std::fill(_states_so2.begin(), _states_so2.end(), Eigen::Rotation2Df(0.0)); };

    /* Construct with an array of angles. assumes radians input. */
    explicit State(const Eigen::Vector4f& angles) {
        std::transform(angles.begin(), angles.end(), _states_so2.begin(),
                       [](const auto angle) { return Eigen::Rotation2Df(angle); });
    }

    /* convert the so2 state vector to an eiger euler angle vector and return */
    Eigen::Vector4f to_angles() const {
        Eigen::Vector4f angles;
        std::transform(_states_so2.begin(), _states_so2.end(), angles.begin(),
                       [](const auto& rot_mat) { return rot_mat.angle(); });
        return angles;
    }

    /* in-place box plus with the perturbation R4 vector */
    void box_plus(const Perturbation& perturbation) {
        std::transform(_states_so2.begin(), _states_so2.end(), perturbation.begin(),
                       _states_so2.begin(),
                       [](const Eigen::Rotation2Df& state, const float pert_angle) {
                           return exponential_map(pert_angle) * state;
                       });
    }
    /* non-const reference to an element of the so2 state vector  */
    Eigen::Rotation2Df& operator[](const std::size_t idx) {
        assert(idx < 4 && "Index must be less than 4");
        return _states_so2.at(idx);
    }
    /* const reference to an element of the so2 state vector */
    const Eigen::Rotation2Df& operator[](const std::size_t idx) const {
        assert(idx < 4 && "Index must be less than 4");
        return _states_so2.at(idx);
    }
    void printStateAngles() const {
        for (const auto& angle : to_angles()) {
            std::cout << "Angle: " << angle << " radians\n";
        }
    }

private:
    std::array<Eigen::Rotation2Df, 4> _states_so2;
};

auto calculate_error_and_jacobian(const State& state, const Measurement& measurement);
};  // namespace rotsync
