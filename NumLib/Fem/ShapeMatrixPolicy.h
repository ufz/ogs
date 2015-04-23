/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SHAPE_MATRIX_POLICY_H_
#define SHAPE_MATRIX_POLICY_H_

#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>

/// An implementation of ShapeMatrixPolicy using fixed size (compile-time) eigen
/// matrices and vectors.
template <typename ShapeFunction>
struct EigenFixedShapeMatrixPolicy
{
    template <int N, int M>
    using _MatrixType = Eigen::Matrix<double, N, M, Eigen::RowMajor>;
    template <int N>
    using _VectorType = Eigen::Matrix<double, N, 1>;

    using NodalMatrixType = _MatrixType<ShapeFunction::NPOINTS, ShapeFunction::NPOINTS>;
    using NodalVectorType = _VectorType<ShapeFunction::NPOINTS>;
    using DimNodalMatrixType = _MatrixType<ShapeFunction::DIM, ShapeFunction::NPOINTS>;
    using DimMatrixType = _MatrixType<ShapeFunction::DIM, ShapeFunction::DIM>;

    using ShapeMatrices =
        NumLib::ShapeMatrices<
            NodalVectorType,
            DimNodalMatrixType,
            DimMatrixType,
            DimNodalMatrixType>;
};

/// An implementation of ShapeMatrixPolicy using dynamic size eigen matrices and
/// vectors.
template <typename ShapeFunction>
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

    using ShapeMatrices =
        NumLib::ShapeMatrices<
            NodalVectorType,
            DimNodalMatrixType,
            DimMatrixType,
            DimNodalMatrixType>;
};

/// Default choice of the ShapeMatrixPolicy.
template <typename ShapeFunction>
using ShapeMatrixPolicyType = EigenFixedShapeMatrixPolicy<ShapeFunction>;
//using ShapeMatrixPolicyType = EigenDynamicShapeMatrixPolicy<ShapeFunction>;

#endif  // OGS_USE_EIGEN

//static_assert(std::is_class<ShapeMatrixPolicyType<>::value,
        //"ShapeMatrixPolicyType was not defined.");

#endif  // SHAPE_MATRIX_POLICY_H_
