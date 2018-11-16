# Look for csparse; note the difference in the directory specifications!
FIND_PATH(LBFGS_INCLUDE_DIR NAMES lbfgs.h
  PATHS
  /usr/include/suitesparse
  /usr/include
  /opt/local/include
  /usr/local/include
  /sw/include
  /usr/include/ufsparse
  /opt/local/include/ufsparse
  /usr/local/include/ufsparse
  /sw/include/ufsparse
  )

FIND_LIBRARY(LBFGS_LIBRARY NAMES lbfgs
  PATHS
  /usr/lib
  /usr/local/lib
  /opt/local/lib
  /sw/lib
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LBFGS DEFAULT_MSG
  LBFGS_INCLUDE_DIR LBFGS_LIBRARY)