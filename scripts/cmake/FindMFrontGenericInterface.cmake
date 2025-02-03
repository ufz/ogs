# - Try to find MGIS
# Once done this will define
#  MFrontGenericInterface_FOUND - System has Mgis
#  MFrontGenericInterface_INCLUDE_DIRS - The Mgis include directories
#  MFrontGenericInterface_LIBRARY_DIRS - The library directories needed to use Mgis
#  MFrontGenericInterface_LIBRARIES    - The libraries needed to use Mgis

if (MGIS_DIR OR MFrontGenericInterface_DIR)
  # in cache already
  SET(MFrontGenericInterface_FIND_QUIETLY TRUE)
endif (MGIS_DIR OR MFrontGenericInterface_DIR)

STRING(REPLACE ":" ";" MYSEARCH_PATH "$ENV{CMAKE_PREFIX_PATH}")
find_file(MFrontGenericInterface_CONFIG_FILE MFrontGenericInterfaceConfig.cmake
  PATHS ${MYSEARCH_PATH}
  PATHS ${MGIS_DIR}/build $ENV{MGIS_DIR}/build ${MGIS_DIR} $ENV{MGIS_DIR} $ENV{MFrontGenericInterface_DIR} ${MFrontGenericInterface_DIR}
  PATH_SUFFIXES share/mgis/cmake
  NO_DEFAULT_PATH
  DOC "The MFrontGenericInterface configuration file")
mark_as_advanced(FORCE MFrontGenericInterface_CONFIG_FILE)


if (MFrontGenericInterface_CONFIG_FILE)
  # Extract the directory name
  get_filename_component(MFrontGenericInterface_CONFIG_PATH ${MFrontGenericInterface_CONFIG_FILE} DIRECTORY)
  set(MFrontGenericInterface_DIR ${MFrontGenericInterface_CONFIG_PATH})
  set(MFrontGenericInterface_DIR ${MFrontGenericInterface_CONFIG_PATH})
  message(STATUS "Found ${MFrontGenericInterface_CONFIG_FILE}")
else(MFrontGenericInterface_CONFIG_FILE)
  message(FATAL_ERROR "MFrontGenericInterface configuration file not found!")
endif(MFrontGenericInterface_CONFIG_FILE)

mark_as_advanced(FORCE MFrontGenericInterface_CONFIG_FILE)
include(${MFrontGenericInterface_CONFIG_FILE})
