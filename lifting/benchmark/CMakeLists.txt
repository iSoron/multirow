add_executable(lifting-benchmark.run src/main.c)
target_link_libraries (lifting-benchmark.run LINK_PUBLIC lifting_static multirow_static)
