set(COMMON_SOURCES
    src/greedy-2d.c
    src/greedy-nd.c
    src/greedy-bsearch.c
    src/greedy.c
    include/infinity/greedy-2d.h
    include/infinity/greedy-nd.h
    include/infinity/greedy-bsearch.h
    include/infinity/greedy.h)

set(TEST_SOURCES
    tests/greedy-2d-test.cpp
    tests/greedy-nd-test.cpp)

add_library(infinity_static ${COMMON_SOURCES})
set_target_properties(infinity_static PROPERTIES OUTPUT_NAME infinity)
target_include_directories (infinity_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(infinity-test.run ${COMMON_SOURCES} ${TEST_SOURCES})
target_link_libraries(infinity-test.run gtest_main multirow_static lifting_static)
