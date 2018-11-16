SET(ISAM_INCLUDE_DIR /usr/local/include/isam)

SET(ISAM_LIB /usr/local/lib)
FIND_LIBRARY(ISAM_LIBRARY
        NAMES libisam.a
        PATHS ${ISAM_LIB})


#include(FindPackageHandleStandardArgs)
#find_package_handle_standard_args(ISAM DEFAULT_MSG ISAM_INCLUDE_DIR ISAM_LIBRARY)
