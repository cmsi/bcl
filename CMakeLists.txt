cmake_minimum_required(VERSION 3.14)
include(cmake/prefix.cmake)
project(exact CXX)
include(cmake/postfix.cmake)
include(cmake/add_iotest.cmake)

add_subdirectory(test)
