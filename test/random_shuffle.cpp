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
#include "bcl/random_shuffle.hpp"

const unsigned int n = 10;
const unsigned int trial1 = 20;
const unsigned int trial2 = 5;

int main() {
  typedef std::vector<unsigned int> vector_type;

  typedef std::mt19937 engine_type;
  engine_type eng(29411);
  std::uniform_real_distribution<> dist;

  std::cout << "generating " << trial1 << " random permutations of "
       << n << " integers [0..." << n - 1 << "]\n";

  for (unsigned int i = 0; i < trial1; ++i) {
    vector_type result(n);
    for (unsigned int j = 0; j < n; ++j) result[j] = j;
    bcl::random_shuffle(result.begin(), result.end(), eng);
    for (unsigned int j = 0; j < n; ++j) {
      std::cout << result[j] << '\t';
    }
    std::cout << std::endl;
  }
}
