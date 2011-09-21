# Source files
SET( HEADERS
	Configure.h.in
	TimeMeasurementBase.h
	RunTimeTimer.h
	CPUTimeTimer.h
)

SET( SOURCES
	RunTimeTimer.cpp
        CPUTimeTimer.cpp
)

# Create the library
ADD_LIBRARY( Base STATIC ${HEADERS} ${SOURCES} )

SET_TARGET_PROPERTIES(Base PROPERTIES LINKER_LANGUAGE CXX)

