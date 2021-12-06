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

#include <cmath> // for std::sqrt

namespace bcl {

class observable {
public:
  observable() : count_(0), sum_(0), esq_(0) {}
  void operator<<(double x) { sum_ += x; esq_ += x * x; ++count_; }
  double mean() const { return (count_ > 0) ? (sum_ / count_) : 0.; }
  double error() const {
    return (count_ > 1) ?
      std::sqrt((esq_ / count_ - mean() * mean()) / (count_ - 1)) : 0.;
  }
private:
  unsigned int count_;
  double sum_, esq_;
};

} // end namespace bcl
