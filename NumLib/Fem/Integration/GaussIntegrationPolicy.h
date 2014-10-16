/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_GAUSSINTEGRATIONPOLICY_H_
#define NUMLIB_GAUSSINTEGRATIONPOLICY_H_

#include "NumLib/Fem/Integration/IntegrationGaussRegular.h"
#include "NumLib/Fem/Integration/IntegrationGaussTri.h"

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
struct GaussIntegrationPolicy<MeshLib::Tri>
{
    using MeshElement = MeshLib::Tri;
    using IntegrationMethod = NumLib::IntegrationGaussTri;
};

}   // namespace NumLib

#endif  // NUMLIB_GAUSSINTEGRATIONPOLICY_H_
