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

static const unsigned int n = 13;
static const unsigned int samples = 10000000;

int main() {
try {
  std::cout << "number of bins = " << n << std::endl;
  std::cout << "number of samples = " << samples << std::endl;

  // random number generator
  typedef std::mt19937 engine_type;
  engine_type eng(29411);
  std::uniform_real_distribution<> dist;

  // generate weights
  std::vector<double> weights(n);
  for (auto& w : weights) w = std::pow(dist(eng), 3.0);
  double total = std::accumulate(weights.begin(), weights.end(), 0.0);
  std::cout << "[weights]\n";
  for (std::size_t i = 0; i < n; ++i) std::cout << weights[i] << ' ';
  std::cout << std::endl;

  // generate markov chain
  int x = 0;
  {
    std::cout << "[Metropolis]\n";
    bcl::markov<engine_type> mc(bcl::metropolis(), weights);
    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t) ++accum[x = mc(x, eng)];
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
    std::cout << "[Metropolis (choose_next)]\n";
    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t)
      ++accum[x = bcl::metropolis::choose_next(weights, x, eng)];
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
    bcl::markov<engine_type> mc(bcl::heatbath(), weights);
    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t) ++accum[x = mc(x, eng)];
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
    std::cout << "[heat bath (choose_next)]\n";
    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t)
      ++accum[x = bcl::heatbath::choose_next(weights, x, eng)];
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
    bcl::markov<engine_type> mc(bcl::st2010(), weights);
    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t) ++accum[x = mc(x, eng)];
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
    std::cout << "[Suwa-Todo 2010 (choose_next)]\n";
    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t)
      ++accum[x = bcl::st2010::choose_next(weights, x, eng)];
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
    std::cout << "[Suwa-Todo 2013]\n";
    bcl::markov<engine_type> mc(bcl::st2013(), weights);
    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t) ++accum[x = mc(x, eng)];
    std::cout << "bin\tweight\t\tresult\t\tdiff\t\tsigma\t\tdiff/sigma\n";
    for (unsigned int i = 0; i < n; ++i) {
      double diff = std::abs((weights[i] / total) - (accum[i] / samples));
      double sigma = std::sqrt(accum[i]) / samples;
      std::cout << i << "\t" << (weights[i] / total) << "    \t"
                << (accum[i] / samples) << "    \t" << diff << "    \t"
                << sigma << "    \t" << (diff / sigma) << std::endl;
    }
  }
}
catch (const std::exception& excp) {
  std::cerr << excp.what() << std::endl;
  std::exit(-1); }
catch (...) {
  std::cerr << "Unknown exception occurred!" << std::endl;
  std::exit(-1); }
}
