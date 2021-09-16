# Find GMP libraries and include path

# Locate libraries and headers for GMP.
find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_path(GMPXX_INCLUDE_DIR NAMES gmpxx.h)
find_library(GMP_LIB NAMES gmp libgmp)
find_library(GMPXX_LIB NAMES gmpxx libgmpxx)

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIB GMPXX_INCLUDE_DIR GMPXX_LIB)

# Determine if the package was found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
  DEFAULT_MSG
  GMP_LIB GMP_INCLUDE_DIR GMPXX_LIB GMPXX_INCLUDE_DIR)
