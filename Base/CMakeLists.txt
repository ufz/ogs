# Source files
SET ( Base_Files
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
	binarySearch.cpp
	DateTools.cpp
        CPUTimeTimer.cpp
	RunTimeTimer.cpp
	StringTools.cpp
)
SOURCE_GROUP( Base FILES ${Base_Files})

# Create the library
ADD_LIBRARY( Base STATIC ${Base_Files})

SET_TARGET_PROPERTIES(Base PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        ../GeoLib
        ../MathLib
        .
)

