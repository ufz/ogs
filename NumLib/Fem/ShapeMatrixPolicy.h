/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SHAPE_MATRIX_POLICY_H_
#define SHAPE_MATRIX_POLICY_H_

#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"

#ifdef OGS_USE_EIGEN
#include <Eigen/Dense>

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
        using type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    };

    template <>
    struct EigenMatrixType<0, 1>
    {
        using type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    };

}   // detail

/// An implementation of ShapeMatrixPolicy using fixed size (compile-time) eigen
/// matrices and vectors.
template <typename ShapeFunction, unsigned GlobalDim>
struct EigenFixedShapeMatrixPolicy
{
    template <int N>
    using _VectorType = typename ::detail::EigenMatrixType<N, 1>::type;

    template <int N, int M>
    using _MatrixType = typename ::detail::EigenMatrixType<N, M>::type;

    using NodalMatrixType = _MatrixType<ShapeFunction::NPOINTS, ShapeFunction::NPOINTS>;
    using NodalVectorType = _VectorType<ShapeFunction::NPOINTS>;
    using DimNodalMatrixType = _MatrixType<ShapeFunction::DIM, ShapeFunction::NPOINTS>;
    using DimMatrixType = _MatrixType<ShapeFunction::DIM, ShapeFunction::DIM>;
    using GlobalDimNodalMatrixType = _MatrixType<GlobalDim, ShapeFunction::NPOINTS>;
    using GlobalDimMatrixType = _MatrixType<GlobalDim, GlobalDim>;

    using ShapeMatrices =
        NumLib::ShapeMatrices<
            NodalVectorType,
            DimNodalMatrixType,
            DimMatrixType,
            GlobalDimNodalMatrixType>;
};

/// An implementation of ShapeMatrixPolicy using dynamic size eigen matrices and
/// vectors.
template <typename ShapeFunction, unsigned GlobalDim>
struct EigenDynamicShapeMatrixPolicy
{
    // Dynamic size local matrices are much slower in allocation than their
    // fixed counterparts.

     using _MatrixType =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
     using _VectorType =
        Eigen::Matrix<double, Eigen::Dynamic, 1>;

    using NodalMatrixType = _MatrixType;
    using NodalVectorType = _VectorType;
    using DimNodalMatrixType = _MatrixType;
    using DimMatrixType = _MatrixType;
    using GlobalDimNodalMatrixType = _MatrixType;
    using GlobalDimMatrixType = _MatrixType;

    using ShapeMatrices =
        NumLib::ShapeMatrices<
            NodalVectorType,
            DimNodalMatrixType,
            DimMatrixType,
            GlobalDimNodalMatrixType>;
};

/// Default choice of the ShapeMatrixPolicy.
template <typename ShapeFunction, unsigned GlobalDim>
using ShapeMatrixPolicyType = EigenFixedShapeMatrixPolicy<ShapeFunction, GlobalDim>;
//using ShapeMatrixPolicyType = EigenDynamicShapeMatrixPolicy<ShapeFunction, GlobalDim>;

#endif  // OGS_USE_EIGEN

//static_assert(std::is_class<ShapeMatrixPolicyType<>::value,
        //"ShapeMatrixPolicyType was not defined.");

#endif  // SHAPE_MATRIX_POLICY_H_
