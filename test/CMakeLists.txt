#  Copyright Synge Todo 2012.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${Boost_ROOT_DIR})

FOREACH(name heatbath markov markov-2 metropolis random_choice random_shuffle st2010 st2013)
  add_executable(${name} ${name}.cpp)
  add_alps_test(${name})
ENDFOREACH(name)
