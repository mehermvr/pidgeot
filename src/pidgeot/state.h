#pragma once
#include <pb_utils/numbers.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <iostream>
#include <ranges>
#include <unordered_map>
#include <utility>

namespace pidgeot {

struct AtomicState {
  // index is not really needed currently. can be removed. only other use is in
  // construction of State from a vector of Atoms
  int index;
  Eigen::Rotation2Dd rotation;
  bool fixed = false;

  AtomicState(const int idx, const double angle) : index(idx), rotation(Eigen::Rotation2Dd(angle)) {}
  AtomicState(const int idx, const Eigen::Rotation2Dd rot) : index(idx), rotation(rot) {}

  /* overloading for right hand side multiplication */
  friend AtomicState operator*(AtomicState lhs, const Eigen::Rotation2Dd& transform) {
    lhs.rotation = lhs.rotation * transform;
    return lhs;
  }
  /* overloading for left hand side multiplication */
  friend AtomicState operator*(const Eigen::Rotation2Dd& transform, AtomicState rhs) {
    rhs.rotation = transform * rhs.rotation;
    return rhs;
  }

  friend std::ostream& operator<<(std::ostream& os, const AtomicState& state_atom) {
    os << state_atom.index << ": " << state_atom.rotation.smallestPositiveAngle() << "\u00B0";
    return os;
  }
};

class State {
private:
  std::unordered_map<int, AtomicState> _elements;

public:
  /* construct states of length a certain size with default angle 0.0 */
  explicit State(int size) {
    std::ranges::for_each(std::views::iota(0, size), [this](const int idx) {
      _elements.emplace(idx, AtomicState{idx, 0.0});
    });
  }

  explicit State(const std::vector<double>& angles) {
    std::ranges::for_each(std::views::enumerate(angles), [this](const auto& enum_elem) {
      const auto& [idx, angle] = enum_elem;
      _elements.emplace(idx, AtomicState{pb_utils::saturate_cast<int>(idx), angle});
    });
  }

  explicit State(const std::vector<AtomicState>& state_atoms) {
    std::ranges::for_each(state_atoms,
                          [this](const auto& state_atom) { _elements.emplace(state_atom.index, state_atom); });
  }

  /* thin wrapper around the map */
  auto begin() { return _elements.begin(); }
  auto begin() const { return _elements.cbegin(); }
  auto end() { return _elements.end(); }
  auto end() const { return _elements.cend(); }
  auto size() const { return _elements.size(); }

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
      auto pert_rotation = Eigen::Rotation2Dd(pert_angle);   // exponential map
      _elements.at(idx) = pert_rotation * _elements.at(idx); // left multiplication
    });
  }

  // fix states to indicate to optimization routines not to move them. doesn't enforce
  // it
  void fix_state(int idx) { _elements.at(idx).fixed = true; }

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
