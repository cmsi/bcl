cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(bcl)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

find_package(BoostSrc)

enable_testing()
include(CTest)
include(add_alps_test)
add_subdirectory(test)
