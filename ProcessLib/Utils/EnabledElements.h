/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/Fem/FiniteElement/ElementTraitsLagrange.h"

#ifndef OGS_MAX_ELEMENT_DIM
static_assert(false, "The macro OGS_MAX_ELEMENT_DIM is undefined.");
#endif

#ifndef OGS_MAX_ELEMENT_ORDER
static_assert(false, "The macro OGS_MAX_ELEMENT_ORDER is undefined.");
#endif

namespace ProcessLib
{
namespace detail
{
static constexpr bool enable_element_type_simplex =
#ifdef OGS_ENABLE_ELEMENT_SIMPLEX
    true;
#else
    false;
#endif

static constexpr bool enable_element_type_cuboid =
#ifdef OGS_ENABLE_ELEMENT_CUBOID
    true;
#else
    false;
#endif

static constexpr bool enable_element_type_prism =
#ifdef OGS_ENABLE_ELEMENT_PRISM
    true;
#else
    false;
#endif

static constexpr bool enable_element_type_pyramid =
#ifdef OGS_ENABLE_ELEMENT_PYRAMID
    true;
#else
    false;
#endif

// Faces of tets, pyramids and prisms are triangles
static constexpr bool enable_element_type_tri = enable_element_type_simplex ||
                                                enable_element_type_pyramid ||
                                                enable_element_type_prism;

// Faces of hexes, pyramids and prisms are quads
static constexpr bool enable_element_type_quad = enable_element_type_cuboid ||
                                                 enable_element_type_pyramid ||
                                                 enable_element_type_prism;

using ZeroOrOneD = std::tuple<MeshLib::Point, MeshLib::Line, MeshLib::Line3>;

using Cuboids = std::tuple<MeshLib::Quad, MeshLib::Quad8, MeshLib::Quad9,
                           MeshLib::Hex, MeshLib::Hex20>;

using Simplices =
    std::tuple<MeshLib::Tri, MeshLib::Tri6, MeshLib::Tet, MeshLib::Tet10>;

using Prisms = std::tuple<MeshLib::Prism, MeshLib::Prism15>;

using Pyramids = std::tuple<MeshLib::Pyramid, MeshLib::Pyramid13>;

using Triangles = std::tuple<MeshLib::Tri, MeshLib::Tri6>;

using Quads = std::tuple<MeshLib::Quad, MeshLib::Quad8, MeshLib::Quad9>;

/**
 * Determines if the given element is contained in the given element group and
 * if furthermore the element group is enabled.
 */
template <typename Elements, typename Element>
constexpr bool isElementEnabledImpl(bool is_group_enabled)
{
    return BaseLib::TMP::contains<Elements, Element>() && is_group_enabled;
}

auto constexpr isElementEnabled = []<typename ElementTraits>(ElementTraits*)
{
    using Element = typename ElementTraits::Element;

    if constexpr (ElementTraits::ShapeFunction::ORDER > OGS_MAX_ELEMENT_ORDER)
    {
        return false;
    }
    if constexpr (Element::dimension > OGS_MAX_ELEMENT_DIM)
    {
        return false;
    }

    return isElementEnabledImpl<ZeroOrOneD, Element>(true) ||
           isElementEnabledImpl<Cuboids, Element>(enable_element_type_cuboid) ||
           isElementEnabledImpl<Simplices, Element>(
               enable_element_type_simplex) ||
           isElementEnabledImpl<Prisms, Element>(enable_element_type_prism) ||
           isElementEnabledImpl<Pyramids, Element>(
               enable_element_type_pyramid) ||
           isElementEnabledImpl<Triangles, Element>(enable_element_type_tri) ||
           isElementEnabledImpl<Quads, Element>(enable_element_type_quad);
};

}  // namespace detail

/// List of all element types for which the local assemblers of OGS processes
/// will be compiled.
using EnabledElementTraitsLagrange =
    decltype(BaseLib::TMP::filter<NumLib::AllElementTraitsLagrange>(
        detail::isElementEnabled));
}  // namespace ProcessLib
