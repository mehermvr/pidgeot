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

namespace pidgeot {

struct AtomicState {
  int index;
  Eigen::Rotation2Dd rotation;

  AtomicState(const int idx, const double angle) : index(idx), rotation(Eigen::Rotation2Dd(angle)) {}
  AtomicState(const int idx, const Eigen::Rotation2Dd rot) : index(idx), rotation(rot) {}

  /* overloading for left hand side multiplication */
  AtomicState& operator*=(const Eigen::Rotation2Dd& transform) {
    this->rotation = transform * this->rotation;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, const AtomicState& state_atom) {
    os << state_atom.index << ": " << lutils::rad2deg(state_atom.rotation.angle()) << "\u00B0";
    return os;
  }
};

class State {
private:
  std::unordered_map<int, AtomicState> _elements;

public:
  /* construct states of length a certain size with default angle 0.0 */
  explicit State(int size) {
    for (int idx : std::views::iota(0, size)) {
      _elements.emplace(std::piecewise_construct, std::forward_as_tuple(idx), std::forward_as_tuple(idx, 0.0));
    }
  }

  explicit State(const std::vector<double>& angles) {
    for (const auto& [idx, angle] : std::views::enumerate(angles)) {
      _elements.emplace(std::piecewise_construct, std::forward_as_tuple(idx), std::forward_as_tuple(idx, angle));
    }
  }

  explicit State(const std::vector<AtomicState>& state_atoms) {
    for (const auto& state_atom : state_atoms) {
      _elements.emplace(state_atom.index, state_atom);
    }
  }
  /* thin wrapper around the map */
  auto begin() { return _elements.begin(); }
  auto begin() const { return _elements.cbegin(); }
  auto end() { return _elements.end(); }
  auto end() const { return _elements.cend(); }

  auto& operator[](const int idx) { return _elements.at(idx); }
  const auto& operator[](const int idx) const { return _elements.at(idx); }

  /* box plus with the perturbation vector */
  /* templating to handle eigen vectors or std vectors */
  template <typename T>
    requires std::ranges::input_range<T>
  void box_plus(const T& perturbation) {
    assert(perturbation.size() == _elements.size() && "Size of perturbation must be same as state");
    std::ranges::for_each(std::views::enumerate(perturbation), [this](const auto& enum_elem) {
      const auto& [idx, pert_angle] = enum_elem;
      auto pert_rotation = Eigen::Rotation2Dd(pert_angle); // exponential map
      /* this is fine */
      /* std::cout << _elements.at(idx); */
      /* this is not fine */
      /* std::cout << _elements[idx]; */
      _elements.at(idx) *= pert_rotation; // left multiplication
    });
  }

  // Overload << operator for printing, expensive since a copy is made for sorting
  friend std::ostream& operator<<(std::ostream& os, const State& state) {
    std::vector<std::pair<int, AtomicState>> pairs(state.begin(), state.end());
    std::ranges::sort(pairs, [](const auto& A, const auto& B) { return A.first < B.first; });
    os << "State {";
    for (const auto& [idx, state_atom] : pairs) {
      os << state_atom << ", ";
    }
    os << "}\n";
    return os;
  }
};

} // namespace pidgeot
