set(COMMON_SOURCES
    src/cg.c
    src/double.c
    src/geometry.c
    src/lfree2d.c
    src/lp.c
    src/mir.c
    src/util.c
    src/rational.c
    src/stats.c
    include/multirow/cg.h
    include/multirow/double.h
    include/multirow/geometry.h
    include/multirow/lfree2d.h
    include/multirow/lp.h
    include/multirow/mir.h
    include/multirow/rational.h
    include/multirow/stats.h
    include/multirow/params.h
    include/multirow/util.h)

set(TEST_SOURCES
    tests/double-test.cpp
    tests/cg-test.cpp
    tests/geometry-test.cpp)

add_library(multirow_static ${COMMON_SOURCES})
set_target_properties(multirow_static PROPERTIES OUTPUT_NAME lifting)
target_link_libraries(multirow_static ${CPLEX_LIBRARIES})
target_include_directories (multirow_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(multirow-test.run ${COMMON_SOURCES} ${TEST_SOURCES})
target_link_libraries(multirow-test.run ${CPLEX_LIBRARIES} gtest_main)
