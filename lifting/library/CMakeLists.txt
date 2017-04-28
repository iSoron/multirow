include_directories(include)

set(COMMON_SOURCES
    src/lifting.c
    src/lifting-mip.c
    include/lifting/lifting-mip.h
    include/lifting/lifting.h)

set(TEST_SOURCES
    tests/lifting-test.cpp)

add_library(lifting_static ${COMMON_SOURCES})
set_target_properties(lifting_static PROPERTIES OUTPUT_NAME lifting)
target_include_directories (lifting_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(lifting-test.run ${COMMON_SOURCES} ${TEST_SOURCES})
target_link_libraries(lifting-test.run gtest_main multirow_static)
