set(COMMON_SOURCES
        src/infinity-2d.c
        src/infinity-nd.c
        src/infinity.c
        include/infinity/infinity-2d.h
        include/infinity/infinity-nd.h
        include/infinity/infinity.h)

set(TEST_SOURCES
        tests/infinity-2d-test.cpp
        tests/infinity-nd-test.cpp
        tests/infinity-test.cpp)

add_library(infinity_static ${COMMON_SOURCES})
set_target_properties(infinity_static PROPERTIES OUTPUT_NAME infinity)
target_include_directories(infinity_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(infinity-test.run ${COMMON_SOURCES} ${TEST_SOURCES})
target_link_libraries(infinity-test.run gtest_main multirow_static lifting_static)
