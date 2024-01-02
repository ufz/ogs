/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Core>

#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"

namespace detail
{
/// Forwards the Eigen::Matrix type for general N and M.
/// There is a partial specialization for M = 1 to store the matrix in
/// column major storage order.
template <int N, int M>
struct EigenMatrixType
{
    using type = Eigen::Matrix<double, N, M, Eigen::RowMajor>;
};

/// Specialization for Nx1 matrices which can be stored only in column major
/// form in Eigen-3.2.5.
template <int N>
struct EigenMatrixType<N, 1>
{
    using type = Eigen::Matrix<double, N, 1, Eigen::ColMajor>;
};

/// Specialization for 0xM matrices. Using fixed size Eigen matrices here
/// would lead to zero sized arrays, which cause compilation errors on
/// some compilers.
template <int M>
struct EigenMatrixType<0, M>
{
    using type =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
};

template <>
struct EigenMatrixType<0, 1>
{
    using type =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
};

}  // namespace detail

/// An implementation of MatrixPolicy using fixed size (compile-time) Eigen
/// matrices and vectors.
struct EigenFixedMatrixPolicy
{
    template <int N>
    using VectorType = typename ::detail::EigenMatrixType<N, 1>::type;

    template <int N>
    using RowVectorType = typename ::detail::EigenMatrixType<1, N>::type;

    template <int N, int M>
    using MatrixType = typename ::detail::EigenMatrixType<N, M>::type;
};

/// An implementation of ShapeMatrixPolicy using fixed size (compile-time) Eigen
/// matrices and vectors.
template <typename ShapeFunction, int GlobalDim>
struct EigenFixedShapeMatrixPolicy : EigenFixedMatrixPolicy
{
    using NodalMatrixType =
        MatrixType<ShapeFunction::NPOINTS, ShapeFunction::NPOINTS>;
    using NodalVectorType = VectorType<ShapeFunction::NPOINTS>;
    using DimVectorType = VectorType<ShapeFunction::DIM>;
    using NodalRowVectorType = RowVectorType<ShapeFunction::NPOINTS>;
    using DimNodalMatrixType =
        MatrixType<ShapeFunction::DIM, ShapeFunction::NPOINTS>;
    using DimMatrixType = MatrixType<ShapeFunction::DIM, ShapeFunction::DIM>;
    using GlobalDimNodalMatrixType =
        MatrixType<GlobalDim, ShapeFunction::NPOINTS>;
    using GlobalDimMatrixType = MatrixType<GlobalDim, GlobalDim>;
    using GlobalDimVectorType = VectorType<GlobalDim>;

    using ShapeMatrices =
        NumLib::ShapeMatrices<NodalRowVectorType, DimNodalMatrixType,
                              DimMatrixType, GlobalDimNodalMatrixType>;
};

/// An implementation of MatrixPolicy using dynamic size Eigen matrices and
/// vectors.
/// \note Dynamic size local matrices are much slower in allocation than their
/// fixed counterparts.
struct EigenDynamicMatrixPolicy
{
    template <int>
    using VectorType = Eigen::Matrix<double, Eigen::Dynamic, 1>;

    template <int>
    using RowVectorType = Eigen::Matrix<double, 1, Eigen::Dynamic>;

    template <int, int>
    using MatrixType =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
};

/// An implementation of ShapeMatrixPolicy using dynamic size Eigen matrices and
/// vectors.
template <typename ShapeFunction, int GlobalDim>
struct EigenDynamicShapeMatrixPolicy : EigenDynamicMatrixPolicy
{
    using NodalMatrixType = MatrixType<0, 0>;
    using NodalVectorType = VectorType<0>;
    using DimVectorType = VectorType<0>;
    using NodalRowVectorType = RowVectorType<0>;
    using DimNodalMatrixType = MatrixType<0, 0>;
    using DimMatrixType = MatrixType<0, 0>;
    using GlobalDimNodalMatrixType = MatrixType<0, 0>;
    using GlobalDimMatrixType = MatrixType<0, 0>;
    using GlobalDimVectorType = VectorType<0>;

    using ShapeMatrices =
        NumLib::ShapeMatrices<NodalRowVectorType, DimNodalMatrixType,
                              DimMatrixType, GlobalDimNodalMatrixType>;
};

#ifdef OGS_EIGEN_DYNAMIC_SHAPE_MATRICES
using MatrixPolicyType = EigenDynamicMatrixPolicy;

template <typename ShapeFunction, int GlobalDim>
using ShapeMatrixPolicyType =
    EigenDynamicShapeMatrixPolicy<ShapeFunction, GlobalDim>;

const unsigned OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG = 1;
#else
using MatrixPolicyType = EigenFixedMatrixPolicy;

template <typename ShapeFunction, int GlobalDim>
using ShapeMatrixPolicyType =
    EigenFixedShapeMatrixPolicy<ShapeFunction, GlobalDim>;

const unsigned OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG = 0;
#endif

// static_assert(std::is_class<ShapeMatrixPolicyType<>::value,
//"ShapeMatrixPolicyType was not defined.");
