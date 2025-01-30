/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on January 21, 2025, 4:10 PM
 */

#pragma once

#include <tuple>

namespace ParameterLib
{
class SpatialPosition;
}
namespace MaterialPropertyLib
{
class Phase;
class VariableArray;

/// It computes fluid density and viscosity for single phase flow model.
std::tuple<double, double> getFluidDensityAndViscosity(
    double const t, double const dt, ParameterLib::SpatialPosition const& pos,
    Phase const& fluid_phase, VariableArray& vars);

/// It computes fluid density for single phase flow model.
double getFluidDensity(double const t, double const dt,
                       ParameterLib::SpatialPosition const& pos,
                       Phase const& fluid_phase, VariableArray& vars);
}  // namespace MaterialPropertyLib
