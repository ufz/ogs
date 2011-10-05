# Source files
SET( HEADERS
	AnalyticalGeometry.h        
	LinearInterpolation.h  
	MathTools.h
	Vector3.h
	EarClippingTriangulation.h  
	LinkedTriangle.h
	max.h
        sparse.h
	vector_io.h
	LinAlg/MatrixBase.h
	LinAlg/VectorNorms.h
	LinAlg/Dense/Matrix.h
        LinAlg/Preconditioner/generateDiagPrecond.h
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
	LinAlg/Sparse/amuxCRS.h
        LinAlg/Sparse/CRSMatrix.h
        LinAlg/Sparse/CRSMatrixPThreads.h
        LinAlg/Sparse/CRSMatrixOpenMP.h
        LinAlg/Sparse/CRSSymMatrix.h
        LinAlg/Sparse/SparseMatrixBase.h
)

SET( SOURCES
	AnalyticalGeometry.cpp        
	LinearInterpolation.cpp  
	MathTools.cpp
	EarClippingTriangulation.cpp  
	LinkedTriangle.cpp
        LinAlg/Sparse/amuxCRS.cpp
        LinAlg/Solvers/BiCGStab.cpp
        LinAlg/Solvers/CG.cpp
        LinAlg/Solvers/GMRes.cpp
	LinAlg/Solvers/GaussAlgorithm.cpp
        LinAlg/Solvers/TriangularSolve.cpp
	LinAlg/Preconditioner/generateDiagPrecond.cpp
)

INCLUDE_DIRECTORIES (
	.
	../Base
	../GeoLib
)

# Create the library
ADD_LIBRARY( MathLib STATIC ${HEADERS} ${SOURCES} )

SET_TARGET_PROPERTIES(MathLib PROPERTIES LINKER_LANGUAGE CXX)

