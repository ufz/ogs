/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Elements/Point.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Pyramid.h"

#include "NumLib/Fem/Integration/IntegrationGaussRegular.h"
#include "NumLib/Fem/Integration/IntegrationGaussTri.h"
#include "NumLib/Fem/Integration/IntegrationGaussTet.h"
#include "NumLib/Fem/Integration/IntegrationGaussPrism.h"
#include "NumLib/Fem/Integration/IntegrationGaussPyramid.h"
#include "NumLib/Fem/Integration/IntegrationPoint.h"

namespace NumLib
{

/// An integration policy providing an integration method suitable for the given
/// mesh element.
/// Gauss-Legendre integration method is used. The default implementation is
/// choosing the regular-element integration, for other elements a
/// specialization must be provided, as for example for triangles
/// (\see GaussIntegrationPolicy<MeshLib::Tri>).
/// The integration method depends on the dimension of the element and correctly
/// chosen number and placement of the integration points within the element.
template <typename MeshElement_>
struct GaussIntegrationPolicy
{
    using MeshElement = MeshElement_;
    using IntegrationMethod =
        NumLib::IntegrationGaussRegular<MeshElement::dimension>;
};

template<>
struct GaussIntegrationPolicy<MeshLib::Point>
{
    using MeshElement = MeshLib::Point;
    using IntegrationMethod = NumLib::IntegrationPoint;
};

template<>
struct GaussIntegrationPolicy<MeshLib::Tri>
{
    using MeshElement = MeshLib::Tri;
    using IntegrationMethod = NumLib::IntegrationGaussTri;
};

template<>
struct GaussIntegrationPolicy<MeshLib::Tri6>
{
    using MeshElement = MeshLib::Tri6;
    using IntegrationMethod = NumLib::IntegrationGaussTri;
};

template<>
struct GaussIntegrationPolicy<MeshLib::Tet>
{
    using MeshElement = MeshLib::Tri;
    using IntegrationMethod = NumLib::IntegrationGaussTet;
};

template<>
struct GaussIntegrationPolicy<MeshLib::Tet10>
{
    using MeshElement = MeshLib::Tet10;
    using IntegrationMethod = NumLib::IntegrationGaussTet;
};

template<>
struct GaussIntegrationPolicy<MeshLib::Prism>
{
    using MeshElement = MeshLib::Prism;
    using IntegrationMethod = NumLib::IntegrationGaussPrism;
};

template<>
struct GaussIntegrationPolicy<MeshLib::Prism15>
{
    using MeshElement = MeshLib::Prism15;
    using IntegrationMethod = NumLib::IntegrationGaussPrism;
};

template<>
struct GaussIntegrationPolicy<MeshLib::Pyramid>
{
    using MeshElement = MeshLib::Pyramid;
    using IntegrationMethod = NumLib::IntegrationGaussPyramid;
};

template<>
struct GaussIntegrationPolicy<MeshLib::Pyramid13>
{
    using MeshElement = MeshLib::Pyramid13;
    using IntegrationMethod = NumLib::IntegrationGaussPyramid;
};

}   // namespace NumLib
