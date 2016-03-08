/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cassert>
#include "UnifiedMatrixSetters.h"

#ifdef OGS_USE_EIGEN

// Dense Eigen matrix/vector //////////////////////////////////////////

namespace MathLib
{

void setMatrix(Eigen::MatrixXd& m,
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

void setMatrix(Eigen::MatrixXd& m, Eigen::MatrixXd const& tmp)
{
    m = tmp;
}

void addToMatrix(Eigen::MatrixXd& m,
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

double norm(Eigen::VectorXd const& x) { return x.norm(); }

void setVector(Eigen::VectorXd& v, std::initializer_list<double> values)
{
    assert((std::size_t) v.size() == values.size());
    auto it = values.begin();
    for (std::size_t i=0; i<values.size(); ++i) v[i] = *(it++);
}

} // namespace MathLib

#endif // OGS_USE_EIGEN


#ifdef USE_PETSC

// Global PETScMatrix/PETScVector //////////////////////////////////////////

#include "MathLib/LinAlg/PETSc/PETScVector.h"

namespace MathLib
{

double norm(PETScVector const& x)
{
    return x.getNorm();
}

void setVector(PETScVector& v,
               std::initializer_list<double> values)
{
    (void) v; (void) values;
    // TODO implement
}


void setMatrix(PETScMatrix& m,
               PETScMatrix::IndexType const rows,
               PETScMatrix::IndexType const cols,
               std::initializer_list<double> values)
{
    (void) m; (void) rows; (void) cols; (void) values;
    // TODO implement
}

void setMatrix(PETScMatrix& m, Eigen::MatrixXd const& tmp)
{

    (void) m; (void) tmp;
    // TODO implement
}

void addToMatrix(PETScMatrix& m,
                 PETScMatrix::IndexType const rows,
                 PETScMatrix::IndexType const cols,
                 std::initializer_list<double> values)
{
    (void) m; (void) rows; (void) cols; (void) values;
    // TODO implement
}

} // namespace MathLib


#elif defined(OGS_USE_EIGEN)

// Sparse global EigenMatrix/EigenVector //////////////////////////////////////////

#include "MathLib/LinAlg/Eigen/EigenVector.h"

namespace MathLib
{

void setVector(EigenVector& v,
                      std::initializer_list<double> values)
{
    setVector(v.getRawVector(), values);
}

void setMatrix(EigenMatrix& m,
               EigenMatrix::IndexType const rows,
               EigenMatrix::IndexType const cols,
               std::initializer_list<double> values)
{
    using IndexType = EigenMatrix::IndexType;
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

void setMatrix(EigenMatrix& m, Eigen::MatrixXd const& tmp)
{
    m.getRawMatrix() = tmp.sparseView();
}

void addToMatrix(EigenMatrix& m,
                 EigenMatrix::IndexType const rows,
                 EigenMatrix::IndexType const cols,
                 std::initializer_list<double> values)
{
    using IndexType = EigenMatrix::IndexType;
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

} // namespace MathLib

#endif // OGS_USE_EIGEN
