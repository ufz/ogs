# Source files
SET( HEADERS
	amuxCRS.h
        CRSMatrix.h
        CRSMatrixPThreads.h
        CRSMatrixOpenMP.h
        CRSSymMatrix.h
        sparse.h
        SparseMatrixBase.h
)

SET( SOURCES
        amuxCRS.cpp
)

INCLUDE_DIRECTORIES (
	.
	${PROJECT_BINARY_DIR}/Base
)

# Create the library
ADD_LIBRARY( MathLib STATIC ${HEADERS} ${SOURCES} )

SET_TARGET_PROPERTIES(MathLib PROPERTIES LINKER_LANGUAGE CXX)

