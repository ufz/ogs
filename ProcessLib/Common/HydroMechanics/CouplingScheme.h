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
    /// scheme. In the HydroMechanics process the fixed-stress split has been
    /// implemented as staggered scheme with a stabilization parameter to be
    /// set. For more details see [user guide -
    /// conventions](https://www.opengeosys.org/docs/userguide/basics/conventions/#fixed-stress-split-for-hydro-mechanical-processes).
    double fixed_stress_stabilization_parameter;

    /// An optional indicator to select whether fixing volume stress rate over
    /// time step or over coupling iteration for the staggered scheme. If it is
    /// not given, volume stress rate is fixed over coupling iteration in the
    /// fixed stress splitting approach in each time step. If set, the volume
    /// stress rate is fixed over time step, e.g.
    /// \f[\dot{\sigma}_v^{n} = \dot{\sigma}_v^{n-1}, \f]
    /// where \f$n\f$ is the time step index.
    /// Otherwise, the volume stress rate is fixed over coupling iteration,
    /// e.g.
    /// \f[\dot{\sigma}_v^{n, k} = \dot{\sigma}_v^{n, k-1}, \f]
    /// where \f$k\f$ is the coupling iteration index, and
    /// \f[ \dot{()}^{n, k} = \left(()^{n, k} - ()^{n-1}\right)/dt. \f]
    bool fixed_stress_over_time_step;
};

using CouplingScheme = std::variant<Monolithic, Staggered>;

CouplingScheme parseCouplingScheme(BaseLib::ConfigTree const& config);
}  // namespace ProcessLib::Common::HydroMechanics
