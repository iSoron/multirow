include_directories(include)

set(COMMON_SOURCES
    src/cplex_helper.cpp
    src/geometry.cpp
    src/gomory_cut_generator.cpp
    src/knapsack2.cpp
    src/mir_cut_generator.cpp
    src/single_row_cut_generator.cpp
    src/stats.cpp
    src/wedge_cut_generator.cpp
    include/onerow/cplex_helper.hpp
    include/onerow/cplex_helper.tpp
    include/onerow/geometry.hpp
    include/onerow/gomory_cut_generator.hpp
    include/onerow/knapsack2.hpp
    include/onerow/mir_cut_generator.hpp
    include/onerow/params.hpp
    include/onerow/single_row_cut_generator.hpp
    include/onerow/stats.hpp
    include/onerow/wedge_cut_generator.hpp)

set(TEST_SOURCES
    tests/knapsack2_test.cpp
    tests/rational_test.cpp
    tests/single_row_generator_test.cpp
    tests/wedge_cut_generator_test.cpp)

add_library(onerow_static ${COMMON_SOURCES})
set_target_properties(onerow_static PROPERTIES OUTPUT_NAME onerow)
target_include_directories (onerow_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(onerow_static qxx_static ${GMP_LIBRARIES} ${CPLEX_LIBRARIES})

add_executable(onerow-test.run ${COMMON_SOURCES} ${TEST_SOURCES})
target_link_libraries(onerow-test.run gtest_main qxx_static ${GMP_LIBRARIES} ${CPLEX_LIBRARIES})
