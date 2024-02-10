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

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "state.h"
namespace rotsync {
inline auto chi_square(const Eigen::Vector4d& error) { return error.squaredNorm(); };

class Solver {
private:
    State _initial_state;
    State _state;
    Measurement _measurement;
    bool _use_analytic_jacobian;
    int _max_iter;
    double _chi_square_thresh;
    int _iter{0};

    Eigen::Vector4d calculate_error(const State& state) {
        auto calc_error_scalar = [](const Eigen::Rotation2Dd& state_1,
                                    const Eigen::Rotation2Dd& state_2, const double relative_rot) {
            Eigen::Matrix2d prediction =
                    state_1.toRotationMatrix().transpose() * state_2.toRotationMatrix();
            Eigen::Matrix2d measurement_so2 = Eigen::Rotation2Dd(relative_rot).toRotationMatrix();
            return (Eigen::Matrix2d::Identity() - prediction.transpose() * measurement_so2).trace();
        };
        Eigen::Vector4d error = Eigen::Vector4d::Zero();
        for (int i = 0; i < error.size(); i++) {
            error[i] = calc_error_scalar(state[i], state[(i + 1) % 4], _measurement[i]);
        }
        return error;
    }
    Eigen::Matrix4d calculate_numerical_jacobian() {
        Eigen::Matrix4d jacobian = Eigen::Matrix4d::Zero();
        auto initial_error = calculate_error(_state);
        double epsilon = 1e-8;
        for (int idx = 0; idx < jacobian.cols(); idx++) {
            Eigen::Vector4d perturbation = Eigen::Vector4d::Zero();
            perturbation[idx] += epsilon;
            auto plussed_state = _state.box_plus(perturbation);
            auto perturbed_error = calculate_error(plussed_state);
            auto delta_error = perturbed_error - initial_error;
            jacobian.col(idx) = delta_error / epsilon;
        }
        return jacobian;
    }
    Eigen::Matrix4d calculate_analytical_jacobian() {
        Eigen::Matrix4d jacobian = Eigen::Matrix4d::Zero();
        auto wrap_index = [](int idx) { return idx % 4; };
        for (int idx = 0; idx < jacobian.cols(); idx++) {
            jacobian(idx, idx) =
                    2 * std::sin(_state[idx].angle() - _state[wrap_index(idx + 1)].angle() +
                                 _measurement[idx]);
            jacobian(idx, wrap_index(idx + 1)) =
                    -2 * std::sin(_state[idx].angle() - _state[wrap_index(idx + 1)].angle() +
                                  _measurement[idx]);
        }
        return jacobian;
    }
    auto calculate_error_and_jacobian() {
        /* calc error */
        Eigen::Vector4d error = calculate_error(_state);
        Jacobian jacobian;
        if (_use_analytic_jacobian) {
            jacobian = calculate_analytical_jacobian();
        } else {
            jacobian = calculate_numerical_jacobian();
        }
        return std::tuple(error, jacobian);
    }

public:
    Solver(const State& initial_state,
           const Measurement& measurement,
           bool use_analytic_jacobian = true,
           int max_iter = 1000,
           double chi_square_thresh = 1e-8)
        : _initial_state(initial_state),
          _state(initial_state),
          _measurement(measurement),
          _use_analytic_jacobian(use_analytic_jacobian),
          _max_iter(max_iter),
          _chi_square_thresh(chi_square_thresh) {}

    auto solve(bool verbose = false) {
        while (_iter < _max_iter) {
            auto [error, jacobian] = calculate_error_and_jacobian();

            Eigen::Matrix4d hessian = jacobian.transpose() * jacobian;
            Eigen::Vector4d rhs = -jacobian.transpose() * error;
            Eigen::Vector4d delta_x = Eigen::Vector4d::Zero();

            /* first measurement is fixed */
            delta_x.tail(3) = hessian.block<3, 3>(1, 1).llt().solve(rhs.tail(3));
            _state = _state.box_plus(delta_x);
            if (verbose) {
                std::cout << "Iter: " << _iter << " and chi_squared = " << chi_square(error)
                          << "\n";
                std::cout << "jacobian is\n"
                          << jacobian << "\nand determinant is " << jacobian.determinant() << " \n";
            }
            if (chi_square(error) < _chi_square_thresh) {
                /* this placed here means one extra iteration though */
                if (verbose) {
                    std::cout << "error is\n" << error << "\n";
                    std::cout << "jacobian is\n"
                              << jacobian << "\nand determinant is " << jacobian.determinant()
                              << " \n";
                    std::cout << "hessian is\n"
                              << hessian << "\nand det is " << hessian.determinant() << " \n";
                }
                break;
            }
            _iter++;
        }
        return _state;
    }
};
};  // namespace rotsync
