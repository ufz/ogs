/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUM_LIB_FEM_FINITE_ELEMENT_LOWERSHAPETABLE_H_
#define NUM_LIB_FEM_FINITE_ELEMENT_LOWERSHAPETABLE_H_

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

namespace NumLib
{
template <typename SF>
struct LowerDim;

template <>
struct LowerDim<NumLib::ShapeLine3>
{
    using type = NumLib::ShapeLine2;
};
template <>
struct LowerDim<NumLib::ShapeQuad8>
{
    using type = NumLib::ShapeQuad4;
};
template <>
struct LowerDim<NumLib::ShapeQuad9>
{
    using type = NumLib::ShapeQuad4;
};
template <>
struct LowerDim<NumLib::ShapeHex20>
{
    using type = NumLib::ShapeHex8;
};
template <>
struct LowerDim<NumLib::ShapeTri6>
{
    using type = NumLib::ShapeTri3;
};
template <>
struct LowerDim<NumLib::ShapeTet10>
{
    using type = NumLib::ShapeTet4;
};
template <>
struct LowerDim<NumLib::ShapePrism15>
{
    using type = NumLib::ShapePrism6;
};
template <>
struct LowerDim<NumLib::ShapePyra13>
{
    using type = NumLib::ShapePyra5;
};

}  // namespace NumLib

#endif  // NUM_LIB_FEM_FINITE_ELEMENT_LOWERSHAPETABLE_H_
