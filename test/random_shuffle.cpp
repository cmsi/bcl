/*****************************************************************************
*
* BCL: Balance Condition Library
*
* Copyright (C) 2009-2013 by Hidemaro Suwa <suwamaro@looper.t.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <bcl/random_shuffle.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <vector>

const unsigned int n = 10;
const unsigned int trial1 = 20;
const unsigned int trial2 = 5;
const boost::uint32_t seed = 23094;

int main() {
  typedef boost::mt19937 rng_type;
  typedef std::vector<unsigned int> vector_type;

  rng_type rng(seed);
  boost::variate_generator<rng_type&, boost::uniform_real<> >
    random(rng, boost::uniform_real<>(0, 1));

  std::cout << "generating " << trial1 << " random permutations of "
       << n << " integers [0..." << n - 1 << "]\n";

  for (unsigned int i = 0; i < trial1; ++i) {
    vector_type result(n);
    for (unsigned int j = 0; j < n; ++j) result[j] = j;
    bcl::random_shuffle(result.begin(), result.end(), random);
    for (unsigned int j = 0; j < n; ++j) {
      std::cout << result[j] << '\t';
    }
    std::cout << std::endl;
  }
}
