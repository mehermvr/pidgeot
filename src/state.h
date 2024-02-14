#pragma once

#include <pb_utils/numbers.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <iostream>

namespace rotsync {
using Perturbation = Eigen::Vector4d;
using Measurement = Eigen::Vector4d;
using Jacobian = Eigen::Matrix4d;

inline Eigen::Rotation2Dd exponential_map(const double angle) { return Eigen::Rotation2Dd(angle); }

class State {
public:
    /* default constructor, all absolute rotations are 0 rad */
    State() { std::fill(_states_so2.begin(), _states_so2.end(), Eigen::Rotation2Dd(0.0)); };

    /* Construct with an array of angles. assumes radians input. */
    template <typename T>
    explicit State(const T& angles) {
        std::transform(angles.cbegin(), angles.cend(), _states_so2.begin(),
                       [](const auto angle) { return Eigen::Rotation2Dd(angle); });
    }

    /* convert the so2 state vector to an eiger euler angle vector and return */
    Eigen::Vector4d to_angles() const {
        Eigen::Vector4d angles;
        std::transform(_states_so2.begin(), _states_so2.end(), angles.begin(),
                       [](const auto& rot_mat) { return rot_mat.angle(); });
        return angles;
    }
    /* does this return a copy or a reference? what happens if i modify? */
    const std::array<Eigen::Rotation2Dd, 4>& as_array() const { return _states_so2; }

    /* box plus with the perturbation R4 vector, return a new state*/
    void box_plus(const Perturbation& perturbation) {
        std::transform(_states_so2.cbegin(), _states_so2.cend(), perturbation.cbegin(),
                       _states_so2.begin(),
                       [](const Eigen::Rotation2Dd& state, const double pert_angle) {
                           return exponential_map(pert_angle) * state;
                       });
    }
    /* non-const reference to an element of the so2 state vector  */
    Eigen::Rotation2Dd& operator[](const std::size_t idx) {
        assert(idx < 4 && "Index must be less than 4");
        return _states_so2.at(idx);
    }
    /* const reference to an element of the so2 state vector */
    const Eigen::Rotation2Dd& operator[](const std::size_t idx) const {
        assert(idx < 4 && "Index must be less than 4");
        return _states_so2.at(idx);
    }
    void printStateAngles() const {
        for (const auto& angle : to_angles()) {
            std::cout << "Angle: " << angle << " radians\n";
        }
    }
    // Overload << operator for printing
    friend std::ostream& operator<<(std::ostream& ostream, const State& state) {
        ostream << "State {";
        const auto& angles = state.to_angles();
        for (const auto& angle : state.to_angles()) {
            ostream << angle * 180 / utils::C_PI << "\u00B0,";
        }
        ostream << "}";
        return ostream;
    }

private:
    std::array<Eigen::Rotation2Dd, 4> _states_so2;
};

}  // namespace rotsync
