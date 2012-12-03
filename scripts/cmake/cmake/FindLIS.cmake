# - Try to find LIS
# Once done, this will define
#
#  LIS_FOUND
#  LIS_INCLUDE_DIRS
#  LIS_LIBRARIES

if (NOT LIS_FOUND)

	include(LibFindMacros)
	
	find_path( LIS_INCLUDE_DIR
		NAMES lis.h
		PATHS ${CMAKE_SOURCE_DIR}/../Libs/precompiled)

	if ( UNIX )
		# Tell if the unix system is on 64-bit base
		if(CMAKE_SIZEOF_VOID_P MATCHES "8")
			find_library(LIS_LIBRARIES
				NAMES lis-64
				PATHS ${CMAKE_SOURCE_DIR}/../Libs/precompiled )	
		else (CMAKE_SIZEOF_VOID_P MATCHES "8")
			find_library(LIS_LIBRARIES
				NAMES lis-32
				PATHS ${CMAKE_SOURCE_DIR}/../Libs/precompiled )	
		endif (CMAKE_SIZEOF_VOID_P MATCHES "8")
	else ( UNIX )			
		find_library(LIS_LIBRARIES
			NAMES lisomp
			PATHS ${CMAKE_SOURCE_DIR}/../Libs/precompiled )	
	endif ( UNIX )

	# Set the include dir variables and the libraries and let libfind_process do the rest.
	# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
	if (NOT LIS_LIBRARIES STREQUAL "LIS_LIBRARIES-NOTFOUND" AND NOT LIS_INCLUDE_DIR STREQUAL "LIS_INCLUDE_DIR-NOTFOUND")
		set(LIS_PROCESS_INCLUDES LIS_INCLUDE_DIR)
		set(LIS_PROCESS_LIBS LIS_LIBRARIES)
		libfind_process(LIS)
	else (NOT LIS_LIBRARIES STREQUAL "LIS_LIBRARIES-NOTFOUND" AND NOT LIS_INCLUDE_DIR STREQUAL "LIS_INCLUDE_DIR-NOTFOUND")
		message (STATUS "Warning: LIS not found!")
	endif (NOT LIS_LIBRARIES STREQUAL "LIS_LIBRARIES-NOTFOUND" AND NOT LIS_INCLUDE_DIR STREQUAL "LIS_INCLUDE_DIR-NOTFOUND")
	
endif (NOT LIS_FOUND)

if(LIS_FOUND)
	INCLUDE_DIRECTORIES( ${LIS_INCLUDE_DIR} )
	message(STATUS "LIS found (include: ${LIS_INCLUDE_DIR})")
endif(LIS_FOUND)
