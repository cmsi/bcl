include_directories(${PROJECT_SOURCE_DIR})

FOREACH(name heatbath markov markov-2 metropolis random_choice random_shuffle st2010 st2013)
  add_executable(${name} ${name}.cpp)
  add_iotest(${name})
ENDFOREACH(name)
