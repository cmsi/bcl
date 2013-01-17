/*****************************************************************************
*
* BCL: Balance Condition Library
*
* Copyright (C) 2011-2012 by Hidemaro Suwa <suwamaro@looper.t.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <bcl.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <vector>

static const unsigned int n = 9;
static const unsigned int samples = 1000000;

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::cout << "number of bins = " << n << std::endl;
  std::cout << "number of samples = " << samples << std::endl;

  // random number generator
  typedef boost::mt19937 engine_type;
  typedef boost::variate_generator<boost::mt19937&, boost::uniform_real<> > generator_type;
  engine_type eng(29411);
  generator_type rng(eng, boost::uniform_real<>());

  // generate weights
  std::vector<double> weights(n);
  std::generate(weights.begin(), weights.end(), rng);
  double total = std::accumulate(weights.begin(), weights.end(), 0.0);
  std::cout << "[weights]\n";
  for (std::size_t i = 0; i < n; ++i) std::cout << weights[i] << ' ';
  std::cout << std::endl;

  // generate markov chain
  int x = 0;
  {
    std::cout << "[Metropolis]\n";
    bcl::markov<generator_type> mc(bcl::metropolis(), weights);
    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t) ++accum[x = mc(x, rng)];
    std::cout << "bin\tweight\t\tresult\t\tdiff\t\tsigma\t\tdiff/sigma\n";
    for (unsigned int i = 0; i < n; ++i) {
      double diff = std::abs((weights[i] / total) - (accum[i] / samples));
      double sigma = std::sqrt(accum[i]) / samples;
      std::cout << i << "\t" << (weights[i] / total) << "    \t"
                << (accum[i] / samples) << "    \t" << diff << "    \t"
                << sigma << "    \t" << (diff / sigma) << std::endl;
    }
  }
  {
    std::cout << "[heat bath]\n";
    bcl::markov<generator_type> mc(bcl::heatbath(), weights);
    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t) ++accum[x = mc(x, rng)];
    std::cout << "bin\tweight\t\tresult\t\tdiff\t\tsigma\t\tdiff/sigma\n";
    for (unsigned int i = 0; i < n; ++i) {
      double diff = std::abs((weights[i] / total) - (accum[i] / samples));
      double sigma = std::sqrt(accum[i]) / samples;
      std::cout << i << "\t" << (weights[i] / total) << "    \t"
                << (accum[i] / samples) << "    \t" << diff << "    \t"
                << sigma << "    \t" << (diff / sigma) << std::endl;
    }
  }
  {
    std::cout << "[Suwa-Todo 2010]\n";
    bcl::markov<generator_type> mc(bcl::st2010(), weights);
    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t) ++accum[x = mc(x, rng)];
    std::cout << "bin\tweight\t\tresult\t\tdiff\t\tsigma\t\tdiff/sigma\n";
    for (unsigned int i = 0; i < n; ++i) {
      double diff = std::abs((weights[i] / total) - (accum[i] / samples));
      double sigma = std::sqrt(accum[i]) / samples;
      std::cout << i << "\t" << (weights[i] / total) << "    \t"
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
