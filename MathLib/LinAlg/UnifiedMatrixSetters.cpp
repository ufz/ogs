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
               std::initializer_list<double> values)
{
    using IndexType = Eigen::MatrixXd::Index;

    auto const rows = m.rows();
    auto const cols = m.cols();

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
                 std::initializer_list<double> values)
{
    using IndexType = Eigen::MatrixXd::Index;

    auto const rows = m.rows();
    auto const cols = m.cols();

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

void setVector(Eigen::VectorXd& v, MatrixVectorTraits<Eigen::VectorXd>::Index const index,
               double const value)
{
    v[index] = value;
}

} // namespace MathLib

#endif // OGS_USE_EIGEN


#ifdef USE_PETSC

// Global PETScMatrix/PETScVector //////////////////////////////////////////

#include <numeric>
#include "MathLib/LinAlg/PETSc/PETScVector.h"
#include "MathLib/LinAlg/PETSc/PETScMatrix.h"

namespace MathLib
{

double norm(PETScVector const& x, MathLib::VecNormType nmtype = MathLib::VecNormType::NORM2)
{
    NormType petsc_norm = NORM_1;
    switch(nmtype)
    {
        case MathLib::VecNormType::NORM1:
            petsc_norm = NORM_1;
            break;
        case MathLib::VecNormType::NORM2:
            petsc_norm = NORM_2;
            break;
        case MathLib::VecNormType::INFINITY_N:
            petsc_norm = NORM_INFINITY;
            break;
        default:
            break;
    }

    PetscScalar norm = 0.;
    VecNorm(x.getRawVector(), petsc_norm, &norm);
    return norm;
}

void setVector(PETScVector& v,
               std::initializer_list<double> values)
{
    std::vector<double> const vals(values);
    std::vector<PETScVector::IndexType> idcs(vals.size());
    std::iota(idcs.begin(), idcs.end(), 0);

    v.set(idcs, vals);
}

void setVector(PETScVector& v, MatrixVectorTraits<PETScVector>::Index const index,
               double const value)
{
    v.set(index, value); // TODO handle negative indices
}

void setMatrix(PETScMatrix& m,
               std::initializer_list<double> values)
{
    m.setZero();
    addToMatrix(m, values);
}

void setMatrix(PETScMatrix& m, Eigen::MatrixXd const& tmp)
{
    using IndexType = PETScMatrix::IndexType;

    auto const rows = tmp.rows();
    auto const cols = tmp.cols();

    assert(rows == m.getNumberOfRows() && cols == m.getNumberOfColumns());

    m.setZero();
    std::vector<IndexType> row_idcs(rows);
    std::vector<IndexType> col_idcs(cols);

    std::iota(row_idcs.begin(), row_idcs.end(), 0);
    std::iota(col_idcs.begin(), col_idcs.end(), 0);

    // PETSc wants row-major
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tmp_ = tmp;

    m.add(row_idcs, col_idcs, tmp_);
}

void addToMatrix(PETScMatrix& m,
                 std::initializer_list<double> values)
{
    using IndexType = PETScMatrix::IndexType;

    auto const rows = m.getNumberOfRows();
    auto const cols = m.getNumberOfColumns();

    assert((IndexType) values.size() == rows*cols);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tmp(rows, cols);

    auto it = values.begin();
    for (IndexType r=0; r<rows; ++r) {
        for (IndexType c=0; c<cols; ++c) {
            tmp(r, c) = *(it++);
        }
    }

    std::vector<IndexType> row_idcs(rows);
    std::vector<IndexType> col_idcs(cols);

    std::iota(row_idcs.begin(), row_idcs.end(), 0);
    std::iota(col_idcs.begin(), col_idcs.end(), 0);

    m.add(row_idcs, col_idcs, tmp);
}

} // namespace MathLib


#elif defined(OGS_USE_EIGEN)

// Sparse global EigenMatrix/EigenVector //////////////////////////////////////////

#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"

namespace MathLib
{

void setVector(EigenVector& v,
                      std::initializer_list<double> values)
{
    setVector(v.getRawVector(), values);
}

void setVector(EigenVector& v, MatrixVectorTraits<EigenVector>::Index const index,
               double const value)
{
    v.getRawVector()[index] = value;
}


void setMatrix(EigenMatrix& m,
               std::initializer_list<double> values)
{
    auto const rows = m.getNumberOfRows();
    auto const cols = m.getNumberOfColumns();

    assert(values.size() == rows*cols);
    Eigen::MatrixXd tmp(rows, cols);

    auto it = values.begin();
    for (std::size_t r=0; r<rows; ++r) {
        for (std::size_t c=0; c<cols; ++c) {
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
                 std::initializer_list<double> values)
{
    auto const rows = m.getNumberOfRows();
    auto const cols = m.getNumberOfColumns();

    assert(values.size() == rows*cols);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tmp(rows, cols);

    auto it = values.begin();
    for (std::size_t r=0; r<rows; ++r) {
        for (std::size_t c=0; c<cols; ++c) {
            tmp(r, c) = *(it++);
        }
    }

    m.getRawMatrix() += tmp.sparseView();
}

} // namespace MathLib

#endif // OGS_USE_EIGEN
