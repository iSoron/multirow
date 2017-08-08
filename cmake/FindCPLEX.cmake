set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH}
  $ENV{HOME}/lib/
  /opt/cplex-12.6/cplex/lib/x86-64_linux/static_pic
  /software/cplex-12.6/distribution/cplex/lib/x86-64_linux/static_pic
  /Users/axavier/Applications/IBM/ILOG/CPLEX_Studio1262/cplex/lib/x86-64_osx/static_pic)

find_library(CPLEX_LIBRARIES
  NAMES cplex cplex1220 cplex1240 cplex1260 cplex1261 cplex1262)

find_path(CPLEX_INCLUDE_DIR NAMES ilcplex/cplex.h PATHS
  $ENV{HOME}/include/
  /opt/cplex-12.6/cplex/include/
  /Users/axavier/Applications/IBM/ILOG/CPLEX_Studio1262/cplex/include/
  /software/cplex-12.6/distribution/cplex/include/)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CPLEX
  DEFAULT_MSG
  CPLEX_LIBRARIES
  CPLEX_INCLUDE_DIR)

mark_as_advanced(
  CPLEX_LIBRARIES
  CPLEX_INCLUDE_DIR)
