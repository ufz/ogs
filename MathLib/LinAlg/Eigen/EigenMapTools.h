/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>
#include <vector>
#include <Eigen/Core>

namespace MathLib
{
/*! Creates an Eigen mapped matrix having its entries stored in the given
 * \c data vector.
 *
 * \return An Eigen mapped matrix of the given size. All values of the matrix
 * are set to zero.
 *
 * \pre The passed \c data vector must have zero size.
 * \post The \c data has size \c rows * \c cols.
 *
 * \note The data vector will have the same storage order (row or column major)
 * as the requested matrix type.
 */
template <typename Matrix>
Eigen::Map<Matrix> createZeroedMatrix(std::vector<double>& data,
                                      Eigen::MatrixXd::Index rows,
                                      Eigen::MatrixXd::Index cols)
{
    static_assert(Matrix::IsRowMajor || Matrix::IsVectorAtCompileTime,
                  "The default storage order in OGS is row major storage for "
                  "dense matrices.");
    assert(Matrix::RowsAtCompileTime == Eigen::Dynamic ||
           Matrix::RowsAtCompileTime == rows);
    assert(Matrix::ColsAtCompileTime == Eigen::Dynamic ||
           Matrix::ColsAtCompileTime == cols);
    assert(data.empty());  // in order that resize fills the vector with zeros.

    data.resize(rows * cols);
    return {data.data(), rows, cols};
}

/*! Creates an Eigen mapped matrix having its entries stored in the given
 * \c data vector.
 *
 * This is a convienence method which makes the specification of dynamically
 * allocated Eigen matrices as return type easier.
 */
inline Eigen::Map<
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
createZeroedMatrix(std::vector<double>& data,
                   Eigen::MatrixXd::Index rows,
                   Eigen::MatrixXd::Index cols)
{
    return createZeroedMatrix<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        data, rows, cols);
}

/*! Creates an Eigen mapped matrix from the given data vector.
 *
 * \attention The data vector must have the same storage order (row or column
 * major) as the requested matrix type.
 */
template <typename Matrix>
Eigen::Map<const Matrix> toMatrix(std::vector<double> const& data,
                                  Eigen::MatrixXd::Index rows,
                                  Eigen::MatrixXd::Index cols)
{
    static_assert(Matrix::IsRowMajor || Matrix::IsVectorAtCompileTime,
                  "The default storage order in OGS is row major storage for "
                  "dense matrices.");
    assert(Matrix::RowsAtCompileTime == Eigen::Dynamic ||
           Matrix::RowsAtCompileTime == rows);
    assert(Matrix::ColsAtCompileTime == Eigen::Dynamic ||
           Matrix::ColsAtCompileTime == cols);
    assert(static_cast<Eigen::MatrixXd::Index>(data.size()) == rows * cols);

    return {data.data(), rows, cols};
}

/*! Creates an Eigen mapped matrix from the given data vector.
 *
 * This is a convienence method which makes the specification of dynamically
 * allocated Eigen matrices as return type easier.
 */
inline Eigen::Map<
    const Eigen::
        Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
toMatrix(std::vector<double> const& data,
         Eigen::MatrixXd::Index rows,
         Eigen::MatrixXd::Index cols)
{
    return toMatrix<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        data, rows, cols);
}

/*! Creates an Eigen mapped matrix from the given data vector.
 *
 * \attention The data vector must have the same storage order (row or column
 * major) as the requested matrix type.
 */
template <typename Matrix>
Eigen::Map<Matrix> toMatrix(std::vector<double>& data,
                            Eigen::MatrixXd::Index rows,
                            Eigen::MatrixXd::Index cols)
{
    static_assert(Matrix::IsRowMajor || Matrix::IsVectorAtCompileTime,
                  "The default storage order in OGS is row major storage for "
                  "dense matrices.");
    assert(Matrix::RowsAtCompileTime == Eigen::Dynamic ||
           Matrix::RowsAtCompileTime == rows);
    assert(Matrix::ColsAtCompileTime == Eigen::Dynamic ||
           Matrix::ColsAtCompileTime == cols);
    assert(static_cast<Eigen::MatrixXd::Index>(data.size()) == rows * cols);

    return {data.data(), rows, cols};
}

/*! Creates an Eigen mapped matrix from the given data vector.
 *
 * This is a convienence method which makes the specification of dynamically
 * allocated Eigen matrices as return type easier.
 */
inline Eigen::Map<
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
toMatrix(std::vector<double>& data,
         Eigen::MatrixXd::Index rows,
         Eigen::MatrixXd::Index cols)
{
    return toMatrix<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        data, rows, cols);
}

/*! Creates an Eigen mapped vector having its entries stored in the given
 * \c data vector.
 *
 * \return An Eigen mapped vector of the given size. All values of the vector
 * are set to zero.
 *
 * \pre The passed \c data vector must have zero size.
 * \post The \c data has size \c size.
 */
template <typename Vector>
Eigen::Map<Vector> createZeroedVector(std::vector<double>& data,
                                      Eigen::VectorXd::Index size)
{
    static_assert(Vector::IsVectorAtCompileTime, "A vector type is required.");
    assert(Vector::SizeAtCompileTime == Eigen::Dynamic ||
           Vector::SizeAtCompileTime == size);
    assert(data.empty());  // in order that resize fills the vector with zeros.

    data.resize(size);
    return {data.data(), size};
}

//! Creates an Eigen mapped vector from the given data vector.
template <typename Vector>
Eigen::Map<const Vector> toVector(std::vector<double> const& data,
                                  Eigen::VectorXd::Index size)
{
    static_assert(Vector::IsVectorAtCompileTime, "A vector type is required.");
    assert(Vector::SizeAtCompileTime == Eigen::Dynamic ||
           Vector::SizeAtCompileTime == size);
    assert(static_cast<Eigen::VectorXd::Index>(data.size()) == size);

    return {data.data(), size};
}

//! Creates an Eigen mapped vector from the given data vector.
template <typename Vector>
Eigen::Map<Vector> toVector(std::vector<double>& data,
                            Eigen::VectorXd::Index size)
{
    static_assert(Vector::IsVectorAtCompileTime, "A vector type is required.");
    assert(Vector::SizeAtCompileTime == Eigen::Dynamic ||
           Vector::SizeAtCompileTime == size);
    assert(static_cast<Eigen::VectorXd::Index>(data.size()) == size);

    return {data.data(), size};
}

/*! Creates an Eigen mapped vector from the given data vector.
 *
 * This is a convienence method which makes the specification of dynamically
 * allocated Eigen vectors as return type easier.
 */
inline Eigen::Map<const Eigen::VectorXd> toVector(
    std::vector<double> const& data)
{
    return {data.data(), static_cast<Eigen::VectorXd::Index>(data.size())};
}

/*! Creates an Eigen mapped vector from the given data vector.
 *
 * This is a convienence method which makes the specification of dynamically
 * allocated Eigen vectors as return type easier.
 */
inline Eigen::Map<Eigen::VectorXd> toVector(
    std::vector<double>& data)
{
    return {data.data(), static_cast<Eigen::VectorXd::Index>(data.size())};
}

}  // MathLib
