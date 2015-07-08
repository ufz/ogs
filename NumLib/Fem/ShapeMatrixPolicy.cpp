/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NumLib/Fem/ShapeMatrixPolicy.h"

#ifdef OGS_USE_EIGEN

template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeHex20  , 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeHex8   , 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapePrism15, 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapePrism6 , 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapePyra13 , 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapePyra5  , 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad4  , 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad8  , 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad9  , 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTet10  , 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTet4   , 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTri3   , 3>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTri6   , 3>;

template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 2>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 2>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad4  , 2>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad8  , 2>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad9  , 2>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTri3   , 2>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeTri6   , 2>;

template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 1>;
template struct EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 1>;

template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeHex20  , 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeHex8   , 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapePrism15, 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapePrism6 , 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapePyra13 , 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapePyra5  , 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad4  , 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad8  , 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad9  , 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTet10  , 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTet4   , 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri3   , 3>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri6   , 3>;

template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 2>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 2>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad4  , 2>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad8  , 2>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad9  , 2>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri3   , 2>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri6   , 2>;

template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 1>;
template struct EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 1>;


#endif  // OGS_USE_EIGEN
