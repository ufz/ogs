/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TemplateIsoparametric.h"

template class NumLib::TemplateIsoparametric<NumLib::ShapeHex20  , EigenFixedShapeMatrixPolicy<NumLib::ShapeHex20  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeHex8   , EigenFixedShapeMatrixPolicy<NumLib::ShapeHex8   , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeLine2  , EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeLine3  , EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapePrism15, EigenFixedShapeMatrixPolicy<NumLib::ShapePrism15, 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapePrism6 , EigenFixedShapeMatrixPolicy<NumLib::ShapePrism6 , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapePyra13 , EigenFixedShapeMatrixPolicy<NumLib::ShapePyra13 , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapePyra5  , EigenFixedShapeMatrixPolicy<NumLib::ShapePyra5  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad4  , EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad4  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad8  , EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad8  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad9  , EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad9  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeTet10  , EigenFixedShapeMatrixPolicy<NumLib::ShapeTet10  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeTet4   , EigenFixedShapeMatrixPolicy<NumLib::ShapeTet4   , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeTri3   , EigenFixedShapeMatrixPolicy<NumLib::ShapeTri3   , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeTri6   , EigenFixedShapeMatrixPolicy<NumLib::ShapeTri6   , 3>>;

template class NumLib::TemplateIsoparametric<NumLib::ShapeLine2  , EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 2>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeLine3  , EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 2>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad4  , EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad4  , 2>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad8  , EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad8  , 2>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad9  , EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad9  , 2>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeTri3   , EigenFixedShapeMatrixPolicy<NumLib::ShapeTri3   , 2>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeTri6   , EigenFixedShapeMatrixPolicy<NumLib::ShapeTri6   , 2>>;

template class NumLib::TemplateIsoparametric<NumLib::ShapeLine2  , EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 1>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeLine3  , EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 1>>;

template class NumLib::TemplateIsoparametric<NumLib::ShapeHex20  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeHex20  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeHex8   , EigenDynamicShapeMatrixPolicy<NumLib::ShapeHex8   , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeLine2  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeLine3  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapePrism15, EigenDynamicShapeMatrixPolicy<NumLib::ShapePrism15, 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapePrism6 , EigenDynamicShapeMatrixPolicy<NumLib::ShapePrism6 , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapePyra13 , EigenDynamicShapeMatrixPolicy<NumLib::ShapePyra13 , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapePyra5  , EigenDynamicShapeMatrixPolicy<NumLib::ShapePyra5  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad4  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad4  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad8  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad8  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad9  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad9  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeTet10  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeTet10  , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeTet4   , EigenDynamicShapeMatrixPolicy<NumLib::ShapeTet4   , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeTri3   , EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri3   , 3>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeTri6   , EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri6   , 3>>;

template class NumLib::TemplateIsoparametric<NumLib::ShapeLine2  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 2>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeLine3  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 2>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad4  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad4  , 2>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad8  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad8  , 2>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad9  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad9  , 2>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeTri3   , EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri3   , 2>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeTri6   , EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri6   , 2>>;

template class NumLib::TemplateIsoparametric<NumLib::ShapeLine2  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 1>>;
template class NumLib::TemplateIsoparametric<NumLib::ShapeLine3  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 1>>;
