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

static const unsigned int n = 8;

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::cout << "number of bins = " << n << std::endl;

  // random number generator
  boost::mt19937 eng(29411);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    rng(eng, boost::uniform_real<>());

  // generate weights
  std::vector<double> weights(n);
  std::generate(weights.begin(), weights.end(), rng);
  std::cout << "[weights]\n";
  for (std::size_t i = 0; i < n; ++i) std::cout << weights[i] << ' ';
  std::cout << std::endl;

  // generate transtition matrix
  std::vector<std::vector<double> > transition_matrix;
  bcl::st2010::generate_transition_matrix_resize(weights, transition_matrix);

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
