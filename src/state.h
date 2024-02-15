#pragma once
#include <pb_utils/numbers.h>

#include "utils.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <iostream>
#include <ranges>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace pigeotto {
using Perturbation = Eigen::Vector4d;
using Jacobian = Eigen::Matrix4d;

inline Eigen::Rotation2Dd exponential_map(const double angle) { return Eigen::Rotation2Dd(angle); }

struct AtomicState {
  int index;
  Eigen::Rotation2Dd element;

  AtomicState(const int idx, const double angle) : element(Eigen::Rotation2Dd(angle)), index(idx) {}

  friend std::ostream& operator<<(std::ostream& os, const AtomicState& state_atom) {
    os << state_atom.index << ": " << state_atom.element.angle() << "\u00B0";
    return os;
  }
};

class State {
public:
  std::unordered_map<int, AtomicState> elements;

  explicit State(int size) {
    for (int idx : std::views::iota(0, size)) {
      elements.emplace(std::piecewise_construct, std::forward_as_tuple(idx), std::forward_as_tuple(idx, 0.0));
    }
  }

  explicit State(const std::vector<double>& angles) {
    for (const auto& [idx, angle] : std::views::enumerate(angles)) {
      elements.emplace(std::piecewise_construct, std::forward_as_tuple(idx), std::forward_as_tuple(idx, angle));
    }
  }
  explicit State(const std::vector<AtomicState>& state_atoms) {
    for (const auto& state_atom : state_atoms) {
      elements.emplace(state_atom.index, state_atom);
    }
  }
  /* convert the so2 state vector to an eiger euler angle vector and return */
  /* Eigen::Vector4d to_angles() const { */
  /*   Eigen::Vector4d angles; */
  /*   std::ranges::transform(_states, angles.begin(), [](const auto& state_atom) { return rot_mat.angle(); }); */
  /*   return angles; */
  /* } */

  /* box plus with the perturbation R4 vector, return a new state*/
  /* void box_plus(const Perturbation& perturbation) { */
  /*   std::ranges::transform( */
  /*       _states, perturbation, _states.begin(), */
  /*       [](const Eigen::Rotation2Dd& state, const double pert_angle) { return exponential_map(pert_angle) * state;
   * }); */
  /* } */
  /* const reference to an element of the so2 state vector */
  /* const Eigen::Rotation2Dd& operator[](const std::size_t idx) const { */
  /*   assert(idx < 4 && "Index must be less than 4"); */
  /*   return _states.at(idx); */
  /* } */
  // Overload << operator for printing, expensive since a copy is made for sorting
  friend std::ostream& operator<<(std::ostream& os, const State& state) {
    std::vector<std::pair<int, AtomicState>> pairs(state.elements.begin(), state.elements.end());
    std::sort(pairs.begin(), pairs.end(), [](const auto& A, const auto& B) { return A.first < B.first; });
    os << "State {";
    for (const auto& [idx, state_atom] : pairs) {
      os << state_atom << ", ";
    }
    os << "}\n";
    return os;
  }
};

} // namespace pigeotto
