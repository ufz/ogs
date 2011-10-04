# Source files
SET( HEADERS
        binarySearch.h
        Configure.h.in
	CPUTimeTimer.h
        DateTools.h
        FileFinder.h
        FileTools.h
        printList.h
        quicksort.h
	RunTimeTimer.h
        StringTools.h
        swap.h
	TimeMeasurementBase.h
        uniqueListInsert.h
        wait.h
)

SET( SOURCES
	binarySearch.cpp
	DateTools.cpp
        CPUTimeTimer.cpp
	RunTimeTimer.cpp
	StringTools.cpp
)

# Create the library
ADD_LIBRARY( Base STATIC ${HEADERS} ${SOURCES} )

SET_TARGET_PROPERTIES(Base PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        ../GEOLib
        ../MathLib
        .
)

