/*
   Copyright (C) 2009-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>,
                              Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#pragma once

#include <cmath>
#include <numeric>

namespace bcl {

template<class MAT>
bool check_probability_conservation(MAT const& transition_matrix, double tolerance = 1.0e-10) {
  std::size_t n = transition_matrix.size();
  for (std::size_t i = 0; i < n; ++i) {
    double p = 0;
    for (std::size_t j = 0; j < n; ++j) p += transition_matrix[i][j];
    if (std::abs(p - 1) > tolerance) return false;
  }
  return true;
}

template<class VEC, class MAT>
bool check_balance_condition(VEC const& weights, MAT const& transition_matrix,
  double tolerance = 1.0e-10) {
  typedef typename VEC::value_type value_type;
  std::size_t n = weights.size();
  for (std::size_t j = 0; j < n; ++j) {
    value_type w = 0;
    for (std::size_t i = 0; i < n; ++i) w += weights[i] * transition_matrix[i][j];
    if (std::abs(w - weights[j]) > tolerance * weights[j]) return false;
  }
  return true;
}

template<class VEC, class MAT>
bool check_detailed_balance(VEC const& weights, MAT const& transition_matrix,
  double tolerance = 1.0e-10) {
  typedef typename VEC::value_type value_type;
  value_type sum = std::accumulate(weights.begin(), weights.end(), value_type(0));
  std::size_t n = weights.size();
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      if (std::abs(weights[i] * transition_matrix[i][j] - weights[j] * transition_matrix[j][i]) >
          tolerance * sum) return false;
  return true;
}

template<class VEC, class MAT>
bool check(VEC const& weights, MAT const& transition_matrix, double tolerance = 1.0e-10) {
  return check_probability_conservagion(transition_matrix, tolerance) &&
    check_balance_condition(weights, transition_matrix, tolerance);
}

template<class VEC, class MAT>
typename VEC::value_type average_rejection(VEC const& weights, MAT const& transition_matrix) {
  typename VEC::value_type vs = 0;
  typename VEC::value_type ws = 0;
  std::size_t n = weights.size();
  for (std::size_t i = 0; i < n; ++i) {
    vs = weights[i] * transition_matrix[i][i];
    ws = weights[i];
  }
  return vs / ws;
}

} // end namespace bcl
