/*
 * \file
 * \brief Cumulative include for all available mesh element types.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Point.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Tri.h"

namespace MeshLib
{

/// A list of types listing all mesh element types supported by OGS.
///
/// This type alias is intended for template metaprogramming et al., not for
/// direct use/instantiation.
using AllElementTypes =
    std::tuple<Point, Line, Line3, Quad, Quad8, Quad9, Hex, Hex20, Tri, Tri6,
               Tet, Tet10, Prism, Prism15, Pyramid, Pyramid13>;

}  // namespace MeshLib
