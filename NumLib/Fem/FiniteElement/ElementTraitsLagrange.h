/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/TMP.h"
#include "MeshLib/Elements/Elements.h"
#include "NumLib/Fem/FiniteElement/LowerDimShapeTable.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#include "NumLib/Fem/ShapeFunction/ShapePoint1.h"
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
namespace detail
{
template <typename ShapeFunction, typename Enabled = void>
struct LowerOrderShapeFunctionOrSame
{
    static_assert(ShapeFunction::ORDER < 2,
                  "Only shape functions of order 1 should use this fallback. "
                  "Order 0 is a special case for 0D elements.");
    using type = ShapeFunction;
};

template <typename ShapeFunction>
struct LowerOrderShapeFunctionOrSame<
    ShapeFunction,
    std::void_t<typename NumLib::LowerDim<ShapeFunction>::type>>
{
    using type = typename NumLib::LowerDim<ShapeFunction>::type;
};

template <typename ShapeFunction_>
struct ShapeFunctionTraits
{
    using ShapeFunction = ShapeFunction_;
    using LowerOrderShapeFunction =
        typename LowerOrderShapeFunctionOrSame<ShapeFunction>::type;
};
}  // namespace detail

template <typename Element>
struct ElementTraitsLagrange;

#define OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(ELEMENT, SHAPE_FUNCTION)    \
    template <>                                                            \
    struct ElementTraitsLagrange<MeshLib::ELEMENT>                         \
        : detail::ShapeFunctionTraits<NumLib::SHAPE_FUNCTION>              \
    {                                                                      \
        using Element = MeshLib::ELEMENT;                                  \
        static_assert(                                                     \
            std::is_same_v<Element, typename ShapeFunction::MeshElement>); \
    }

// points and lines
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Point, ShapePoint1);
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Line, ShapeLine2);
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Line3, ShapeLine3);
// quads and hexahedra
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Quad, ShapeQuad4);
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Quad8, ShapeQuad8);
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Quad9, ShapeQuad9);
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Hex, ShapeHex8);
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Hex20, ShapeHex20);
// simplices
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Tri, ShapeTri3);
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Tri6, ShapeTri6);
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Tet, ShapeTet4);
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Tet10, ShapeTet10);
// prisms
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Prism, ShapePrism6);
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Prism15, ShapePrism15);
// pyramids
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Pyramid, ShapePyra5);
OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE(Pyramid13, ShapePyra13);

#undef OGS_SPECIALIZE_ELEMENT_TRAITS_LAGRANGE

using AllElementTraitsLagrange =
    BaseLib::TMP::Map_t<ElementTraitsLagrange, MeshLib::AllElementTypes>;
}  // namespace NumLib
