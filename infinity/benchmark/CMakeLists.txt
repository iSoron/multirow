add_executable(infinity-benchmark.run src/main.c)
target_link_libraries (infinity-benchmark.run LINK_PUBLIC
    lifting_static
    infinity_static
    multirow_static
    m
    pthread)
