find_library(LAPACKE_LIBRARIES
        NAMES lapacke)

find_path(LAPACKE_INCLUDE_DIR NAMES lapacke.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACKE
        DEFAULT_MSG
        LAPACKE_LIBRARIES
        LAPACKE_INCLUDE_DIR)
mark_as_advanced(LAPACKE_LIBRARIES LAPACKE_INCLUDE_DIR)
