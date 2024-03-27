/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Core>
#include <Eigen/LU>
#include <type_traits>

#include "BaseLib/Error.h"

namespace MathLib
{
/// Namespace for vectorized second-order tensor types and operations on those.
/// The tensors are for example, deformation gradients.
namespace VectorizedTensor
{
/// Vectorized tensor size for given displacement dimension.
constexpr int size(int const displacement_dim)
{
    if (displacement_dim == 1)
    {
        return 3;
    }
    if (displacement_dim == 2)
    {
        return 5;
    }
    if (displacement_dim == 3)
    {
        return 9;
    }
    OGS_FATAL(
        "Cannot convert displacement dimension {} to vectorized tensor size.",
        displacement_dim);
}

/// Displacement dimension of a vectorized tensor.
template <typename VectorizedTensor>
constexpr int dimension()
{
    static_assert(VectorizedTensor::ColsAtCompileTime == 1);
    constexpr int rows = VectorizedTensor::RowsAtCompileTime;
    if (rows == 3)
    {
        return 1;
    }
    if (rows == 5)
    {
        return 2;
    }
    if (rows == 9)
    {
        return 3;
    }
    OGS_FATAL(
        "Cannot convert vectorized tensor of size {} to displacement "
        "dimension.",
        rows);
}

/// Vectorized tensor type for given displacement dimension.
/// \note The Eigen vector is always a fixed size vector in contrast to a shape
/// matrix policy types like GMatrixPolicyType::GradientVectorType.
template <int DisplacementDim>
using Type = Eigen::Matrix<double, size(DisplacementDim), 1, Eigen::ColMajor>;

/// Vectorized identity tensor expression. Corresponds to 3x3 identity matrix in
/// all dimensions.
template <int DisplacementDim>
constexpr auto identity()
{
    static_assert(
        0 < DisplacementDim && DisplacementDim <= 3,
        "Identity is implemented only for displacement dimension 1, 2, or 3.");
    return Type<DisplacementDim>::NullaryExpr(
        size(DisplacementDim), 1,
        [](Eigen::Index const row, [[maybe_unused]] Eigen::Index const col)
        {
            assert(col == 0);
            if constexpr (DisplacementDim == 1)
            {  // All entries are diagonal entries.
                return 1;
            }
            if constexpr (DisplacementDim == 2)
            {
                if (row == 0 || row == 3 || row == 4)
                {
                    return 1;
                }
            }
            if constexpr (DisplacementDim == 3)
            {
                if (row == 0 || row == 4 || row == 8)
                {
                    return 1;
                }
            }
            return 0;
        });
}

/// Computes determinant of a vectorized tensor.
template <typename Derived>
double determinant(Eigen::MatrixBase<Derived> const& tensor)
{
    constexpr int displacement_dim = dimension<Derived>();
    static_assert(0 < displacement_dim && displacement_dim <= 3,
                  "Vectorized tensor determinant is implemented only for "
                  "displacement dimension 1, 2, or 3.");

    if constexpr (displacement_dim == 1)
    {
        return tensor[0] * tensor[1] * tensor[2];
    }
    if constexpr (displacement_dim == 2)
    {
        Eigen::Map<Eigen::Matrix2d const> const top_left{
            tensor.derived().data()};

        return top_left.determinant() * tensor[4];
    }
    if constexpr (displacement_dim == 3)
    {
        return Eigen::Map<Eigen::Matrix3d const>(tensor.derived().data())
            .determinant();
    }
}

/// Only a diagonal tensor can be converted to a 1d vectorized tensor.
bool isTensorConvertibleTo1d(Eigen::Matrix3d const& tensor);

/// Only a tensor of form
/// \f$ \left( \begin{array}{ccc}
/// a & b & 0 \\%
/// c & d & 0 \\%
/// 0 & 0 & e \\%
/// \end{array} \right)\f$
/// can be converted to a 2d vectorized tensor.
bool isTensorConvertibleTo2d(Eigen::Matrix3d const& tensor);

/// Converts a 3x3 matrix expression to a vectorized tensor if the conversion
/// for that dimension is possible; see isTensorConvertibleTo1d() and
/// isTensorConvertibleTo2d() functions.
template <int DisplacementDim, typename Derived>
Type<DisplacementDim> toVector(Eigen::MatrixBase<Derived> const& tensor)
{
    static_assert(0 < DisplacementDim && DisplacementDim <= 3,
                  "Conversion to displacement dimension other than 1, 2, or 3 "
                  "is not valid.");

    constexpr int rows = Derived::RowsAtCompileTime;
    constexpr int cols = Derived::ColsAtCompileTime;
    static_assert(rows == 3 || rows == Eigen::Dynamic);
    static_assert(cols == 3 || cols == Eigen::Dynamic);
    if (tensor.rows() != 3 || tensor.cols() != 3)
    {
        OGS_FATAL(
            "Incorrect tensor size, must be 3x3, but tensor is {:d}x{:d}.",
            tensor.rows(), tensor.cols());
    }

    if constexpr (DisplacementDim == 1)
    {
        if (!isTensorConvertibleTo1d(tensor))
        {
            OGS_FATAL(
                "Cannot convert a tensor with non-zero off-diagonal elements "
                "to a 1d vectorized tensor representation.");
        }
        return tensor.diagonal();
    }
    if constexpr (DisplacementDim == 2)
    {
        if (!isTensorConvertibleTo2d(tensor))
        {
            OGS_FATAL(
                "Cannot convert a tensor with non-zero elements at (0, 2), (1, "
                "2), (2, 0), and (2, 1) positions to a 2d vectorized tensor "
                "representation.");
        }
        Type<2> result;
        result.template head<4>() =
            tensor.template block<2, 2>(0, 0).reshaped();
        result(4) = tensor(2, 2);
        return result;
    }
    if constexpr (DisplacementDim == 3)
    {
        return tensor.reshaped();
    }
    OGS_FATAL(
        "Not all cases handled in the VectorizedTensor::toVector() function.");
}

/// Converts a vectorized tensor to a 3x3 matrix.
template <int DisplacementDim>
Eigen::Matrix3d toTensor(Type<DisplacementDim> const& tensor)
{
    static_assert(
        DisplacementDim == 1 || DisplacementDim == 2 || DisplacementDim == 3,
        "Conversion from displacement dimension other than 1, 2, or 3 "
        "is not valid.");

    using Matrix = Eigen::Matrix3d;
    if constexpr (DisplacementDim == 1)
    {
        return tensor.asDiagonal();
    }
    if constexpr (DisplacementDim == 2)
    {
        Matrix m = Matrix::Zero();
        m.template block<2, 2>(0, 0) =
            Eigen::Map<Eigen::Matrix<double, 2, 2> const>(tensor.data());
        m(2, 2) = tensor(4);
        return m;
    }
    if constexpr (DisplacementDim == 3)
    {
        return Eigen::Map<Matrix const>(tensor.data());
    }
}

}  // namespace VectorizedTensor
}  // namespace MathLib
