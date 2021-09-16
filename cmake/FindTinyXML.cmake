# Find TINYXML libraries and include path

# Locate libraries and headers for TinyXML.
find_path(TINYXML_INCLUDE_DIR NAMES tinyxml.h)
find_library(TINYXML_LIB NAMES tinyxml)

mark_as_advanced(TINYXML_INCLUDE_DIR TINYXML_LIB)

# Determine if the package was found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TINYXML
  DEFAULT_MSG
  TINYXML_LIB TINYXML_INCLUDE_DIR)
