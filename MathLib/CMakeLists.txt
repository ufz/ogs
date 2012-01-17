#Source files grouped by a directory
SET ( MathLib_Files
	AnalyticalGeometry.h        
	LinearInterpolation.h  
	MathTools.h
	Vector3.h
	EarClippingTriangulation.h  
	max.h
        sparse.h
	vector_io.h
	AnalyticalGeometry.cpp        
	LinearInterpolation.cpp  
	MathTools.cpp
	EarClippingTriangulation.cpp  
)
SOURCE_GROUP( MathLib FILES ${MathLib_Files})
SET ( SOURCES ${SOURCES} ${MathLib_Files})

SET ( MathLib_LinAlg_Files
	LinAlg/MatrixBase.h
	LinAlg/VectorNorms.h
)
SOURCE_GROUP( MathLib\\LinAlg FILES ${MathLib_LinAlg_Files})
SET ( SOURCES ${SOURCES} ${MathLib_LinAlg_Files})

SET ( MathLib_LinAlg_Dense_Files
	LinAlg/Dense/Matrix.h
)
SOURCE_GROUP( MathLib\\LinAlg\\Dense FILES ${MathLib_LinAlg_Dense_Files})
SET ( SOURCES ${SOURCES} )

SET ( MathLib_LinAlg_Sparse_Files
	LinAlg/Sparse/amuxCRS.h
        LinAlg/Sparse/CRSMatrix.h
        LinAlg/Sparse/CRSMatrixPThreads.h
        LinAlg/Sparse/CRSMatrixOpenMP.h
        LinAlg/Sparse/CRSSymMatrix.h
        LinAlg/Sparse/SparseMatrixBase.h
        LinAlg/Sparse/amuxCRS.cpp
)
SOURCE_GROUP( MathLib\\LinAlg\\Sparse FILES ${MathLib_LinAlg_Sparse_Files})
SET ( SOURCES ${SOURCES} ${MathLib_LinAlg_Sparse_Files})

SET ( MathLib_LinAlg_Solvers_Files
        LinAlg/Solvers/LinearSolver.h
        LinAlg/Solvers/DirectLinearSolver.h
        LinAlg/Solvers/DenseDirectLinearSolver.h
        LinAlg/Solvers/GaussAlgorithm.h
        LinAlg/Solvers/TriangularSolve.h
        LinAlg/Solvers/IterativeLinearSolver.h
        LinAlg/Solvers/solver.h
        LinAlg/Solvers/BiCGStab.h
        LinAlg/Solvers/CG.h
        LinAlg/Solvers/GMRes.h
        LinAlg/Solvers/BiCGStab.cpp
        LinAlg/Solvers/CG.cpp
        LinAlg/Solvers/CGParallel.cpp
        LinAlg/Solvers/GMRes.cpp
	LinAlg/Solvers/GaussAlgorithm.cpp
        LinAlg/Solvers/TriangularSolve.cpp
)
SOURCE_GROUP( MathLib\\LinAlg\\Solvers FILES ${MathLib_LinAlg_Solvers_Files})
SET ( SOURCES ${SOURCES} ${MathLib_LinAlg_Solvers_Files})

SET ( MathLib_LinAlg_Preconditioner_Files
        LinAlg/Preconditioner/generateDiagPrecond.h
	LinAlg/Preconditioner/generateDiagPrecond.cpp
)
SOURCE_GROUP( MathLib\\LinAlg\\Preconditioner FILES ${MathLib_LinAlg_Preconditioner_Files})
SET ( SOURCES ${SOURCES} ${MathLib_LinAlg_Preconditioner_Files})

IF (METIS_FOUND)
	FILE(GLOB MathLib_LinAlg_Sparse_NestedDissection_HEADERS
		LinAlg/Sparse/NestedDissectionPermutation *.h)

	FILE(GLOB MathLib_LinAlg_Sparse_NestedDissection_SOURCES
		LinAlg/Sparse/NestedDissectionPermutation *.cpp)

	SOURCE_GROUP( MathLib\\LinAlg\\Sparse\\NestedDissectionPermutation FILES
		${MathLib_LinAlg_Sparse_NestedDissectionPermutation_HEADERS}
		${MathLib_LinAlg_Sparse_NestedDissectionPermutation_SOURCES}
	)
	SET (SOURCES ${SOURCES} 
		${MathLib_LinAlg_Sparse_NestedDissectionPermutation_HEADERS}
		${MathLib_LinAlg_Sparse_NestedDissectionPermutation_SOURCES}
	)
ENDIF ()

INCLUDE_DIRECTORIES (
	.
	../Base
	../GeoLib
)

IF(METIS_FOUND)
	INCLUDE_DIRECTORIES (
		${METIS_INCLUDE_DIR}
	)
ENDIF()

# Create the library
ADD_LIBRARY( MathLib STATIC 
	${SOURCES}
)

SET_TARGET_PROPERTIES(MathLib PROPERTIES LINKER_LANGUAGE CXX)

