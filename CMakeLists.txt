cmake_minimum_required(VERSION 3.14)
include(cmake/prefix.cmake)
project(exact CXX)
include(cmake/postfix.cmake)

option(BCL_INSTALL "Install bcl." ON)
option(BCL_BUILD_TESTS "Build all bcl's tests." ON)
option(BCL_BUILD_EXAMPLES "Build all bcl's examples." ON)

include(cmake/add_iotest.cmake)

if(BCL_BUILD_TESTS)
  add_subdirectory(test)
endif(BCL_BUILD_TESTS)
