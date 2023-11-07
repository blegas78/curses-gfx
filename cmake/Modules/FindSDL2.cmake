# This module defines
#  SDL2_INCLUDE_DIR, where to find SDL2 header files location, etc.
#  SDL2_FOUND, If false, do not try to use SDL2.
# also defined, but not for general use are
#  SDL2_LIBRARY, where to find the SDL2 library.

FIND_PATH(
  SDL2_INCLUDE_DIR SDL2
  PATHS
    ${SDL2_ROOT}
  PATH_SUFFIXES 
    include
    include/SDL2
  NO_DEFAULT_PATH )

FIND_PATH(
  SDL2_INCLUDE_DIR SDL2
  PATHS
    ${SDL2_ROOT}
  PATH_SUFFIXES 
    include
    include/SDL2 )

SET(SDL2_NAMES ${SDL2_NAMES} SDL2)

FIND_LIBRARY(
  SDL2_LIBRARY 
  NAMES 
    ${SDL2_NAMES}
  PATHS
    ${SDL2_ROOT}
  PATH_SUFFIXES
    /lib
  NO_DEFAULT_PATHS )

FIND_LIBRARY(
  SDL2_LIBRARY 
  NAMES 
    ${SDL2_NAMES}
  PATHS
    ${SDL2_ROOT}
  PATH_SUFFIXES
    /lib )

# handle the QUIETLY and REQUIRED arguments and set SDL2_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
#FIND_PACKAGE_HANDLE_STANDARD_ARGS(SDL2 DEFAULT_MSG SDL2_LIBRARY SDL2_INCLUDE_DIR)
find_package_handle_standard_args(SDL2  DEFAULT_MSG  SDL2_LIBRARY SDL2_INCLUDE_DIR)
#if(SDL2_FOUND)
#	add_definitions(-DSDL2_FOUND)
#endif()
