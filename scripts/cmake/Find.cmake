######################
### Find tools     ###
######################

# Find dot tool from graphviz
FIND_PROGRAM(DOT_TOOL_PATH dot DOC "Dot tool from graphviz")

# Find doxygen
FIND_PACKAGE(Doxygen)

# Find gnu profiler gprof
FIND_PROGRAM(GPROF_PATH gprof DOC "GNU profiler gprof")

FIND_PACKAGE(cppcheck)

######################
### Find libraries ###
######################

# Clang does not have OpenMP support atm, see https://github.com/ufz/ogs/issues/8
IF(NOT COMPILER_IS_CLANG)
	FIND_PACKAGE(OpenMP)
ENDIF () # !clang
IF(OPENMP_FOUND)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

FIND_PACKAGE(Metis)

## Qt4 library ##
IF(NOT OGS_DONT_USE_QT)
	FIND_PACKAGE( Qt4 4.5)
ENDIF(NOT OGS_DONT_USE_QT)

IF ( QT4_FOUND )
	# Enable more modules
	SET(QT_USE_QTOPENGL TRUE)
	SET(QT_USE_QTSQL TRUE)
	SET(QT_USE_QTTEST TRUE)
	SET(QT_USE_QTXML TRUE)
	IF(QT_QTXMLPATTERNS_FOUND)
		SET(QT_USE_QTXMLPATTERNS TRUE)
	ENDIF(QT_QTXMLPATTERNS_FOUND)
	INCLUDE( ${QT_USE_FILE} )
	ADD_DEFINITIONS(${QT_DEFINITIONS})
ENDIF (QT4_FOUND )

## pthread ##
SET ( CMAKE_THREAD_PREFER_PTHREAD ON )
FIND_PACKAGE ( Threads )
IF ( CMAKE_USE_PTHREADS_INIT )
	SET (HAVE_PTHREADS TRUE)
ENDIF (CMAKE_USE_PTHREADS_INIT )

# blas
FIND_PACKAGE ( BLAS )

# lapack
FIND_PACKAGE ( LAPACK )