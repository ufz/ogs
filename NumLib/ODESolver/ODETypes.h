#pragma once

#include <initializer_list>
#include <cassert>

#define USE_EIGEN_PLAIN



#ifdef USE_EIGEN_PLAIN
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

using IndexType = int;
using ODEMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, IndexType>;
using ODEVector = Eigen::VectorXd;

inline void oneShotLinearSolve(ODEMatrix& A, ODEVector& rhs, ODEVector& x)
{
    Eigen::SparseLU<ODEMatrix> slv;
    slv.compute(A);
    x = slv.solve(rhs);
}

inline double norm(ODEVector const& x) { return x.norm(); }

#else
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenLinearSolver.h"

using IndexType = int;
using ODEMatrix = MathLib::EigenMatrix;
using ODEVector = MathLib::EigenVector;

inline void oneShotLinearSolve(ODEMatrix& A, ODEVector& rhs, ODEVector& x)
{
    MathLib::EigenLinearSolver slv(A);
    slv.solve(rhs, x);
}

#endif



inline void setMatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& m,
                      IndexType const rows, IndexType const cols,
                      std::initializer_list<double> values)
{
    assert((IndexType) values.size() == rows*cols);
    Eigen::MatrixXd tmp(rows, cols);

    auto it = values.begin();
    for (IndexType r=0; r<rows; ++r) {
        for (IndexType c=0; c<cols; ++c) {
            tmp(r, c) = *(it++);
        }
    }

    m = tmp.sparseView();
}

inline void setMatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& m,
                      Eigen::MatrixXd const& tmp)
{
    m = tmp.sparseView();
}

inline void addToMatrix(Eigen::SparseMatrix<double, Eigen::RowMajor>& m,
                        IndexType const rows, IndexType const cols,
                        std::initializer_list<double> values)
{
    assert((IndexType) values.size() == rows*cols);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tmp(rows, cols);

    auto it = values.begin();
    for (IndexType r=0; r<rows; ++r) {
        for (IndexType c=0; c<cols; ++c) {
            tmp(r, c) = *(it++);
        }
    }

    m += tmp.sparseView();
}



#ifndef USE_EIGEN_PLAIN
inline void setMatrix(ODEMatrix& m, IndexType const rows, IndexType const cols,
                      std::initializer_list<double> values)
{
    setMatrix(m.getRawMatrix(), rows, cols, values);
}

inline void setMatrix(ODEMatrix& m, Eigen::MatrixXd const& tmp)
{
    m.getRawMatrix() = tmp.sparseView();
}

inline void addToMatrix(ODEMatrix& m, IndexType const rows, IndexType const cols,
                        std::initializer_list<double> values)
{
    addToMatrix(m.getRawMatrix(), rows, cols, values);
}

#endif


enum class NonlinearSolverTag : bool { Picard, Newton };

enum class ODESystemTag : char
{
    FirstOrderImplicitQuasilinear
};
