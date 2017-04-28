cmake_minimum_required(VERSION 2.8)
project(multirow)

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /Users/axavier/Applications/IBM/ILOG/CPLEX_Studio1262/cplex/lib/x86-64_osx/static_pic)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(CPLEX REQUIRED)
include_directories(${CPLEX_INCLUDE_DIR})

include_directories(multirow/include)
include_directories(lifting/library/include)
include_directories(infinity/library/include)

add_subdirectory(googletest/googletest)
include_directories(${gtest_SOURCE_DIR}/include)

add_subdirectory(multirow)
add_subdirectory(lifting/library)
add_subdirectory(lifting/benchmark)
add_subdirectory(infinity/library)
add_subdirectory(infinity/benchmark)
