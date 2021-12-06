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

// ST2013: Rejection optimized reversible kernel (symmetrized version of ST2010)

#pragma once

#include <vector>
#include <algorithm>

namespace bcl {

class st2013 {
public:
  // Generate transition matrix
  template<class VEC, class MAT>
  static void generate_transition_matrix(VEC const& weights, MAT& tm) {
    using std::abs;
    std::size_t n = weights.size();
    std::vector<double> accum(n+1, 0);
    double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    double shift = *(std::max_element(weights.begin(), weights.end())) / sum;
    accum[0] = 0;
    for (std::size_t i = 0; i < n; ++i) accum[i+1] = accum[i] + weights[i] / sum;

    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        tm[i][j] = (std::max(std::min(accum[i+1] + shift, accum[j+1]) -
                             std::max(accum[i] + shift, accum[j]), 0.0) +
                    std::max(std::min(accum[i+1] + shift, accum[j+1] + 1) -
                             std::max(accum[i] + shift, accum[j] + 1), 0.0));
      }
    }
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = i; j < n; ++j) {
        double t = (tm[i][j] + tm[j][i]) / 2;
        tm[i][j] = t / (weights[i] / sum);
        tm[j][i] = t / (weights[j] / sum);
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
    
  // template<class VEC, class ENGINE>
  // static std::size_t choose_next(VEC const& weights, std::size_t present, ENGINE& eng){
  //   std::uniform_real_distribution<> dist;
  //   double sum = *(std::max_element(weights.begin(), weights.end()));
  //   sum -= weights[present] * dist(eng);
  //   for (int i = 0; i < weights.size(); ++i) {
  //     int j = (present + i + 1) % weights.size();
  //     if (sum <= weights[j]) return j;
  //     sum -= weights[j];
  //   }
  //   return present;
  // }
};

} // end namespace bcl
