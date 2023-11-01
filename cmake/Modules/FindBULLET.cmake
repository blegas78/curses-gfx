# This module defines
#  BULLET_INCLUDE_DIR, where to find bullet header files location, etc.
#  BULLET_FOUND, If false, do not try to use BULLET.
# also defined, but not for general use are
#  BULLET_LIBRARY, where to find the BULLET library.

FIND_PATH(
  BULLET_INCLUDE_DIR bullet
  PATHS
    ${BULLET_ROOT}
  PATH_SUFFIXES 
    include
    include/bullet
  NO_DEFAULT_PATH )

FIND_PATH(
  BULLET_INCLUDE_DIR bullet
  PATHS
    ${BULLET_ROOT}
  PATH_SUFFIXES 
    include
    include/bullet )

SET(BULLET_NAMES ${BULLET_NAMES} Bullet Bullet3Collision)

FIND_LIBRARY(
  BULLET_LIBRARY 
  NAMES 
    ${BULLET_NAMES}
  PATHS
    ${BULLET_ROOT}
  PATH_SUFFIXES
    /lib
  NO_DEFAULT_PATHS )

FIND_LIBRARY(
  BULLET_LIBRARY 
  NAMES 
    ${BULLET_NAMES}
  PATHS
    ${BULLET_ROOT}
  PATH_SUFFIXES
    /lib )

# handle the QUIETLY and REQUIRED arguments and set BULLET_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
#FIND_PACKAGE_HANDLE_STANDARD_ARGS(BULLET DEFAULT_MSG BULLET_LIBRARY BULLET_INCLUDE_DIR)
find_package_handle_standard_args(BULLET  DEFAULT_MSG  BULLET_LIBRARY BULLET_INCLUDE_DIR)
#if(BULLET_FOUND)
#	add_definitions(-DBULLET_FOUND)
#endif()
