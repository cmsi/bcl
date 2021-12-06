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

#include <iostream>
#include <random>
#include <vector>
#include "bcl.hpp"

static const unsigned int n = 8;

int main() {
try {
  std::cout << "number of bins = " << n << std::endl;

  // random number generator
  typedef std::mt19937 engine_type;
  engine_type eng(29411);
  std::uniform_real_distribution<> dist;

  // generate weights
  std::vector<double> weights(n);
  for (auto& w : weights) w = dist(eng);
  std::cout << "[weights]\n";
  for (std::size_t i = 0; i < n; ++i) std::cout << weights[i] << ' ';
  std::cout << std::endl;

  // generate transtition matrix
  std::vector<std::vector<double> > transition_matrix;
  bcl::metropolis::generate_transition_matrix_resize(weights, transition_matrix);

  std::cout << "[transition matrix]\n";
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      std::cout << transition_matrix[i][j] << ' ';
    }
    std::cout << std::endl;
  }

  // check transition matrix
  std::cout << "[check transition matrix]\n";
  std::cout << "probability conservation = "
            << (bcl::check_probability_conservation(transition_matrix) ? "pass" : "fail")
            << std::endl;
  std::cout << "detailed balance condition = "
            << (bcl::check_detailed_balance(weights, transition_matrix) ? "pass" : "fail")
            << std::endl;
  std::cout << "balance condition = "
            << (bcl::check_balance_condition(weights, transition_matrix) ? "pass" : "fail")
            << std::endl;

  std::cout << "average rejection = " << bcl::average_rejection(weights, transition_matrix) 
            << std::endl;
}
catch (const std::exception& excp) {
  std::cerr << excp.what() << std::endl;
  std::exit(-1); }
catch (...) {
  std::cerr << "Unknown exception occurred!" << std::endl;
  std::exit(-1); }
}
