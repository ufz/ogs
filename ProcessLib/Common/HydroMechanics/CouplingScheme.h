// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <variant>

namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib::Common::HydroMechanics
{
struct Monolithic
{
};

struct Staggered
{
    /// An optional input to set an algorithmic parameter of the staggered
    /// scheme. In the HydroMechanics and the ThermoHydroMechanics processes
    /// the fixed-stress split has been implemented as staggered scheme with a
    /// stabilization parameter to be set. For more details see [user guide -
    /// conventions](https://www.opengeosys.org/docs/userguide/basics/conventions/#fixed-stress-split-for-hydro-mechanical-processes).
    double fixed_stress_stabilization_parameter;
};

using CouplingScheme = std::variant<Monolithic, Staggered>;

CouplingScheme parseCouplingScheme(BaseLib::ConfigTree const& config);
}  // namespace ProcessLib::Common::HydroMechanics
