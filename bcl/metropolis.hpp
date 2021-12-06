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

#include <algorithm>
#include <cmath>
#include <vector>

namespace bcl {

class metropolis {
public:
  // Generate transition matrix
  template<class VEC, class MAT>
  static void generate_transition_matrix(VEC const& weights, MAT& tm) {
    using std::abs;
    typedef typename VEC::value_type value_type;
    std::size_t n = weights.size();
    for (std::size_t i = 0; i < n; ++i) {
      tm[i][i] = 1;
      if (weights[i] > 0) {
        for (std::size_t j = 0; j < n; ++j) {
          if (i != j) {
            tm[i][j] = std::min(weights[i], weights[j]) / weights[i] / (n-1);
            tm[i][i] -= tm[i][j];
          }
        }
        tm[i][i] = abs(tm[i][i]);
      } else {
        for (std::size_t j = 0; j < n; ++j)
          if (i != j) tm[i][j] = value_type(0);
      }
    }
  }
  template<class VEC, class MAT>
  static void generate_transition_matrix_resize(VEC const& weights, MAT& tm) {
    std::size_t n = weights.size();
    tm.resize(n);
    for (std::size_t i = 0; i < n; ++i) tm[i].resize(n);
    generate_transition_matrix(weights, tm);
  }

  template<class VEC, class ENGINE>
  static std::size_t choose_next(VEC const& weights, std::size_t present, ENGINE& eng) {
    std::uniform_real_distribution<> dist;
    std::size_t proposal = weights.size() * dist(eng);
    return (weights[present] * dist(eng) < weights[proposal]) ? proposal : present;
  }
};

} // end namespace bcl
