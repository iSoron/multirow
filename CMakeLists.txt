cmake_minimum_required(VERSION 2.8)
project(multirow)

set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(CPLEX REQUIRED)
include_directories(${CPLEX_INCLUDE_DIR})

find_package(GMP REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(${gtest_SOURCE_DIR}/include)
include_directories(infinity/library/include)
include_directories(lifting/library/include)
include_directories(multirow/include)
include_directories(onerow/library/include)
include_directories(qxx/include)

add_subdirectory(googletest/googletest)
add_subdirectory(infinity/benchmark)
add_subdirectory(infinity/library)
add_subdirectory(lifting/benchmark)
add_subdirectory(lifting/library)
add_subdirectory(multirow)
add_subdirectory(onerow/benchmark)
add_subdirectory(onerow/library)
add_subdirectory(qxx)
