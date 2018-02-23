/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Elements/Point.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Tri.h"

#include "NumLib/Fem/Integration/IntegrationGaussLegendrePrism.h"
#include "NumLib/Fem/Integration/IntegrationGaussLegendrePyramid.h"
#include "NumLib/Fem/Integration/IntegrationGaussLegendreRegular.h"
#include "NumLib/Fem/Integration/IntegrationGaussLegendreTet.h"
#include "NumLib/Fem/Integration/IntegrationGaussLegendreTri.h"
#include "NumLib/Fem/Integration/IntegrationPoint.h"

namespace NumLib
{
/// An integration policy providing an integration method suitable for the given
/// mesh element.
/// Gauss-Legendre integration method is used. The default implementation is
/// choosing the regular-element integration, for other elements a
/// specialization must be provided, as for example for triangles
/// (\see GaussLegendreIntegrationPolicy<MeshLib::Tri>).
/// The integration method depends on the dimension of the element and correctly
/// chosen number and placement of the integration points within the element.
template <typename MeshElement_>
struct GaussLegendreIntegrationPolicy
{
    using MeshElement = MeshElement_;
    using IntegrationMethod =
        NumLib::IntegrationGaussLegendreRegular<MeshElement::dimension>;
};

template <>
struct GaussLegendreIntegrationPolicy<MeshLib::Point>
{
    using MeshElement = MeshLib::Point;
    using IntegrationMethod = NumLib::IntegrationPoint;
};

template <>
struct GaussLegendreIntegrationPolicy<MeshLib::Tri>
{
    using MeshElement = MeshLib::Tri;
    using IntegrationMethod = NumLib::IntegrationGaussLegendreTri;
};

template <>
struct GaussLegendreIntegrationPolicy<MeshLib::Tri6>
{
    using MeshElement = MeshLib::Tri6;
    using IntegrationMethod = NumLib::IntegrationGaussLegendreTri;
};

template <>
struct GaussLegendreIntegrationPolicy<MeshLib::Tet>
{
    using MeshElement = MeshLib::Tri;
    using IntegrationMethod = NumLib::IntegrationGaussLegendreTet;
};

template <>
struct GaussLegendreIntegrationPolicy<MeshLib::Tet10>
{
    using MeshElement = MeshLib::Tet10;
    using IntegrationMethod = NumLib::IntegrationGaussLegendreTet;
};

template <>
struct GaussLegendreIntegrationPolicy<MeshLib::Prism>
{
    using MeshElement = MeshLib::Prism;
    using IntegrationMethod = NumLib::IntegrationGaussLegendrePrism;
};

template <>
struct GaussLegendreIntegrationPolicy<MeshLib::Prism15>
{
    using MeshElement = MeshLib::Prism15;
    using IntegrationMethod = NumLib::IntegrationGaussLegendrePrism;
};

template <>
struct GaussLegendreIntegrationPolicy<MeshLib::Pyramid>
{
    using MeshElement = MeshLib::Pyramid;
    using IntegrationMethod = NumLib::IntegrationGaussLegendrePyramid;
};

template <>
struct GaussLegendreIntegrationPolicy<MeshLib::Pyramid13>
{
    using MeshElement = MeshLib::Pyramid13;
    using IntegrationMethod = NumLib::IntegrationGaussLegendrePyramid;
};

}  // namespace NumLib
