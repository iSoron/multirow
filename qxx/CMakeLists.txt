set(COMMON_SOURCES
    src/dlu.cpp
    src/dmat.cpp
    src/dvec.cpp
    src/rational.cpp
    src/slu.cpp
    src/smat.cpp
    src/strprintf.cpp
    src/svec.cpp
    include/qxx/dlu.hpp
    include/qxx/dmat.hpp
    include/qxx/dvec.hpp
    include/qxx/rational.hpp
    include/qxx/slu.hpp
    include/qxx/smat.hpp
    include/qxx/strprintf.hpp
    include/qxx/svec.hpp)

add_library(qxx_static ${COMMON_SOURCES})
set_target_properties(qxx_static PROPERTIES OUTPUT_NAME qxx)
target_include_directories (qxx_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
