# Find ARGP libraries and include path

# Locate libraries and headers for ARGP.
find_path(ARGP_INCLUDE_DIR NAMES argp.h)
find_library(ARGP_LIB NAMES argp)

mark_as_advanced(ARGP_INCLUDE_DIR ARGP_LIB)

# Determine if the package was found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARGP
  DEFAULT_MSG
  ARGP_LIB ARGP_INCLUDE_DIR)
