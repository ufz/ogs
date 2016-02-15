/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_ODETYPES_H
#define NUMLIB_ODETYPES_H

// TODO: move that file somewhere else

#include <initializer_list>
#include <cassert>


// Dense Eigen matrix/vector //////////////////////////////////////////
// always enabled

#include <Eigen/LU>


inline void setMatrix(Eigen::MatrixXd& m,
                      Eigen::MatrixXd::Index const rows, Eigen::MatrixXd::Index const cols,
                      std::initializer_list<double> values)
{
    using IndexType = Eigen::MatrixXd::Index;
    assert((IndexType) values.size() == rows*cols);

    auto it = values.begin();
    for (IndexType r=0; r<rows; ++r) {
        for (IndexType c=0; c<cols; ++c) {
            m(r, c) = *(it++);
        }
    }
}

inline void setMatrix(Eigen::MatrixXd& m, Eigen::MatrixXd const& tmp)
{
    m = tmp;
}

inline void addToMatrix(Eigen::MatrixXd& m,
                        Eigen::MatrixXd::Index const rows, Eigen::MatrixXd::Index const cols,
                        std::initializer_list<double> values)
{
    using IndexType = Eigen::MatrixXd::Index;
    assert((IndexType) values.size() == rows*cols);

    auto it = values.begin();
    for (IndexType r=0; r<rows; ++r) {
        for (IndexType c=0; c<cols; ++c) {
            m(r, c) += *(it++);
        }
    }
}

namespace NumLib
{
inline void oneShotLinearSolve(Eigen::MatrixXd& A, Eigen::VectorXd& rhs, Eigen::VectorXd& x)
{
    Eigen::FullPivLU<Eigen::MatrixXd> slv;
    slv.compute(A);
    x = slv.solve(rhs);
}

inline double norm(Eigen::VectorXd const& x) { return x.norm(); }

}

inline void setVector(Eigen::VectorXd& v, std::initializer_list<double> values)
{
    assert((std::size_t) v.size() == values.size());
    auto it = values.begin();
    for (std::size_t i=0; i<values.size(); ++i) v[i] = *(it++);
}


#ifdef USE_PETSC

// Global PETScMatrix/PETScVector //////////////////////////////////////////

#include "MathLib/LinAlg/PETSc/PETScMatrix.h"
#include "MathLib/LinAlg/PETSc/PETScVector.h"


namespace NumLib
{
inline void oneShotLinearSolve(MathLib::PETScMatrix& A, MathLib::PETScVector& rhs,
                               MathLib::PETScVector& x)
{
    (void) A; (void) rhs; (void) x;
    // TODO implement
}

inline double norm(MathLib::PETScVector const& x)
{
    return x.getNorm();
}
}


inline void setVector(MathLib::PETScVector& v,
                      std::initializer_list<double> values)
{
    (void) v; (void) values;
    // TODO implement
}


inline void setMatrix(MathLib::PETScMatrix& m,
                      MathLib::PETScMatrix::IndexType const rows,
                      MathLib::PETScMatrix::IndexType const cols,
                      std::initializer_list<double> values)
{
    (void) m; (void) rows; (void) cols; (void) values;
    // TODO implement
}

inline void setMatrix(MathLib::PETScMatrix& m, Eigen::MatrixXd const& tmp)
{

    (void) m; (void) tmp;
    // TODO implement
}

inline void addToMatrix(MathLib::PETScMatrix& m,
                        MathLib::PETScMatrix::IndexType const rows,
                        MathLib::PETScMatrix::IndexType const cols,
                        std::initializer_list<double> values)
{
    (void) m; (void) rows; (void) cols; (void) values;
    // TODO implement
}


#elif defined(OGS_USE_EIGEN)

// Sparse global EigenMatrix/EigenVector //////////////////////////////////////////

#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenLinearSolver.h"


namespace NumLib
{
inline void oneShotLinearSolve(MathLib::EigenMatrix& A, MathLib::EigenVector& rhs, MathLib::EigenVector& x)
{
    MathLib::EigenLinearSolver slv(A);
    slv.solve(rhs, x);
}
}


inline void setVector(MathLib::EigenVector& v,
                      std::initializer_list<double> values)
{
    setVector(v.getRawVector(), values);
}


inline void setMatrix(MathLib::EigenMatrix& m,
                      MathLib::EigenMatrix::IndexType const rows,
                      MathLib::EigenMatrix::IndexType const cols,
                      std::initializer_list<double> values)
{
    using IndexType = MathLib::EigenMatrix::IndexType;
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

inline void setMatrix(MathLib::EigenMatrix& m, Eigen::MatrixXd const& tmp)
{
    m.getRawMatrix() = tmp.sparseView();
}

inline void addToMatrix(MathLib::EigenMatrix& m,
                        MathLib::EigenMatrix::IndexType const rows,
                        MathLib::EigenMatrix::IndexType const cols,
                        std::initializer_list<double> values)
{
    using IndexType = MathLib::EigenMatrix::IndexType;
    assert((IndexType) values.size() == rows*cols);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tmp(rows, cols);

    auto it = values.begin();
    for (IndexType r=0; r<rows; ++r) {
        for (IndexType c=0; c<cols; ++c) {
            tmp(r, c) = *(it++);
        }
    }

    m.getRawMatrix() += tmp.sparseView();
}

#endif // OGS_USE_EIGEN

#endif // NUMLIB_ODETYPES_H
