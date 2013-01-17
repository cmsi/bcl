/* test/random_choice.C test program
 *
 * Copyright (C) 2012 by Synge Todo <wistaria@comp-phys.org>
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 */

#include <bcl/random_choice.hpp>
#include <boost/random.hpp>
#include <cmath>
#include <vector>

static const unsigned int n = 9;
static const unsigned int samples = 100000;

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
  std::cout << "number of bins = " << n << std::endl;
  std::cout << "number of samples = " << samples << std::endl;

  std::vector<double> weights(n);

  // random number generator
  typedef boost::mt19937 engine_type;
  typedef boost::variate_generator<boost::mt19937&, boost::uniform_real<> > generator_type;
  engine_type eng(29411);
  generator_type rng(eng, boost::uniform_real<>());

  // generate weights
  std::generate(weights.begin(), weights.end(), rng);
  double tw = 0;
  for (unsigned int i = 0; i < n; ++i) tw += weights[i];

  // double-base version
  {
    // random_choice
    bcl::random_choice<generator_type> dist(weights);

    // check
    if (dist.check(weights)) {
      std::cout << "check succeeded\n";
    } else {
      std::cout << "check failed\n";
      std::exit(-1);
    }

    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t) ++accum[dist(rng)];

    std::cout << "bin\tweight\t\tresult\t\tdiff\t\tsigma\t\tdiff/sigma\n";
    for (unsigned int i = 0; i < n; ++i) {
      double diff = std::abs((weights[i] / tw) - (accum[i] / samples));
      double sigma = std::sqrt(accum[i]) / samples;
      std::cout << i << "\t" << (weights[i] / tw) << "    \t"
                << (accum[i] / samples) << "    \t" << diff << "    \t"
                << sigma << "    \t" << (diff / sigma) << std::endl;
    }
  }

  // integer-base version
  {
    // random_choice
    bcl::random_choice<engine_type> dist(weights);
    // bcl::random_choice_walker_i<> dist(weights);

    // check
    if (dist.check(weights)) {
      std::cout << "check succeeded\n";
    } else {
      std::cout << "check failed\n";
      std::exit(-1);
    }

    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t) ++accum[dist(eng)];

    std::cout << "bin\tweight\t\tresult\t\tdiff\t\tsigma\t\tdiff/sigma\n";
    for (unsigned int i = 0; i < n; ++i) {
      double diff = std::abs((weights[i] / tw) - (accum[i] / samples));
      double sigma = std::sqrt(accum[i]) / samples;
      std::cout << i << "\t" << (weights[i] / tw) << "    \t"
                << (accum[i] / samples) << "    \t" << diff << "    \t"
                << sigma << "    \t" << (diff / sigma) << std::endl;
    }
  }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (const std::exception& excp) {
  std::cerr << excp.what() << std::endl;
  std::exit(-1); }
catch (...) {
  std::cerr << "Unknown exception occurred!" << std::endl;
  std::exit(-1); }
#endif
}
