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

#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism15.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism6.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra13.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra5.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad8.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad9.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet10.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet4.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri6.h"

#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>

/// An implementation of ShapeMatrixPolicy using fixed size (compile-time) eigen
/// matrices and vectors.
template <typename ShapeFunction, unsigned GlobalDim>
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

extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeHex20  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeHex8   , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapePrism15, 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapePrism6 , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapePyra13 , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapePyra5  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad4  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad8  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad9  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTet10  , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTet4   , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTri3   , 3>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTri6   , 3>;

extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 2>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 2>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad4  , 2>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad8  , 2>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad9  , 2>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTri3   , 2>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTri6   , 2>;

extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 1>;
extern template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 1>;

extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeHex20  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeHex8   , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapePrism15, 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapePrism6 , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapePyra13 , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapePyra5  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad4  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad8  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad9  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTet10  , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTet4   , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri3   , 3>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri6   , 3>;

extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 2>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 2>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad4  , 2>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad8  , 2>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad9  , 2>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri3   , 2>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri6   , 2>;

extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 1>;
extern template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 1>;
#endif  // OGS_USE_EIGEN

//static_assert(std::is_class<ShapeMatrixPolicyType<>::value,
        //"ShapeMatrixPolicyType was not defined.");

#endif  // SHAPE_MATRIX_POLICY_H_
