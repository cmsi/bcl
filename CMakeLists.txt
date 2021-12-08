cmake_minimum_required(VERSION 3.14)
include(cmake/prefix.cmake)
project(bcl CXX)
include(cmake/postfix.cmake)

option(BCL_INSTALL "Install bcl." ON)
option(BCL_BUILD_TESTS "Build all bcl's tests." ON)
option(BCL_BUILD_EXAMPLES "Build all bcl's examples." ON)

include_directories(${PROJECT_SOURCE_DIR})
include(cmake/add_iotest.cmake)

include(googletest)
FetchContent_MakeAvailable(${FetchContents})

if(BCL_BUILD_TESTS)
  add_subdirectory(test)
endif(BCL_BUILD_TESTS)
