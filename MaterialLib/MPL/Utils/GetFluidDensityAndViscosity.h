// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
