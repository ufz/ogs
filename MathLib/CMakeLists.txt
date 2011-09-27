# Source files
SET( HEADERS
	LinAlg/MatrixBase.h
	LinAlg/Sparse/amuxCRS.h
        LinAlg/Sparse/CRSMatrix.h
        LinAlg/Sparse/CRSMatrixPThreads.h
        LinAlg/Sparse/CRSMatrixOpenMP.h
        LinAlg/Sparse/CRSSymMatrix.h
        sparse.h
        LinAlg/Sparse/SparseMatrixBase.h
)

SET( SOURCES
        LinAlg/Sparse/amuxCRS.cpp
)

INCLUDE_DIRECTORIES (
	.
	${PROJECT_BINARY_DIR}/Base
)

# Create the library
ADD_LIBRARY( MathLib STATIC ${HEADERS} ${SOURCES} )

SET_TARGET_PROPERTIES(MathLib PROPERTIES LINKER_LANGUAGE CXX)

