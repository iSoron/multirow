find_library(CPLEX_LIBRARIES
  NAMES cplex cplex1220 cplex1240 cplex1260 cplex1261 cplex1262)

find_path(CPLEX_INCLUDE_DIR
  NAMES ilcplex/cplex.h
  PATHS /Users/axavier/Applications/IBM/ILOG/CPLEX_Studio1262/cplex/include)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CPLEX
  DEFAULT_MSG
  CPLEX_LIBRARIES
  CPLEX_INCLUDE_DIR)

mark_as_advanced(
  CPLEX_LIBRARIES
  CPLEX_INCLUDE_DIR)
