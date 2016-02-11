#pragma once

#include <initializer_list>
#include <cassert>

#if 0
#include <Eigen/SparseCore>

using IndexType = int;
using Matrix = Eigen::SparseMatrix<double, 0, IndexType>;
using Vector = Eigen::VectorXd;
#else
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenLinearSolver.h"

using IndexType = int;
using Matrix = MathLib::EigenMatrix;
using Vector = MathLib::EigenVector;

using LinearSolver = MathLib::EigenLinearSolver;


inline void setMatrix(Matrix& m, IndexType const rows, IndexType const cols,
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

    m.getRawMatrix() = tmp.sparseView();
}

inline void setMatrix(Matrix& m, Eigen::MatrixXd const& tmp)
{
    m.getRawMatrix() = tmp.sparseView();
}

#endif


enum class NonlinearSolverTag : bool { Picard, Newton };

enum class ODESystemTag : char
{
    FirstOrderImplicitQuasilinear
};
