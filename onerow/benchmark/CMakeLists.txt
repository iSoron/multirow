add_executable(onerow-benchmark.run src/main.cpp)
target_link_libraries (onerow-benchmark.run LINK_PUBLIC
    qxx_static
    onerow_static
    m
    pthread
    ${GMP_LIBRARIES}
    ${CPLEX_LIBRARIES}
    ${OpenMP_LIBRARIES})
