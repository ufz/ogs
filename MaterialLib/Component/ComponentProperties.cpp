/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   ComponentProperties.cpp
 */

#include "ComponentProperties.h"

namespace MaterialLib
{
namespace Component
{
double ComponentProperties::getMolecularDiffusion(
    double const t, ProcessLib::SpatialPosition const& x) const
{
    return _molecular_diffusion(t, x)[0];
}

double ComponentProperties::getLongitudinalDispersivity(
    double const t, ProcessLib::SpatialPosition const& x) const
{
    return _solute_dispersivity_longitudinal(t, x)[0];
}

double ComponentProperties::getTransversalDispersivity(
    double const t, ProcessLib::SpatialPosition const& x) const
{
    return _solute_dispersivity_transverse(t, x)[0];
}
}  // namespace Component
}  // namespace MaterialLib
