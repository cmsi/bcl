cmake_minimum_required(VERSION 3.14)
include(cmake/prefix.cmake)
project(bcl CXX)
include(cmake/postfix.cmake)

include_directories(${PROJECT_SOURCE_DIR})

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(cmake/add_iotest.cmake)
  include(googletest)
endif()
FetchContent_MakeAvailable(${FetchContents})

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  add_subdirectory(test)
endif()
