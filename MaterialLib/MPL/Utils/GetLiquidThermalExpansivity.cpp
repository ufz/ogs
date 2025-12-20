// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "GetLiquidThermalExpansivity.h"

#include "MaterialLib/MPL/Phase.h"

namespace MaterialPropertyLib
{
class Phase;

double getLiquidThermalExpansivity(Phase const& phase,
                                   VariableArray const& vars,
                                   const double density,
                                   ParameterLib::SpatialPosition const& pos,
                                   double const t, double const dt)
{
    // The thermal expansivity is explicitly given in the project file.
    if (phase.hasProperty(
            MaterialPropertyLib::PropertyType::thermal_expansivity))
    {
        return phase
            .property(MaterialPropertyLib::PropertyType::thermal_expansivity)
            .template value<double>(vars, pos, t, dt);
    }

    // The thermal expansivity calculated by the density model directly.
    return (density == 0.0)
               ? 0.0
               : -phase.property(MaterialPropertyLib::PropertyType::density)
                         .template dValue<double>(
                             vars, MaterialPropertyLib::Variable::temperature,
                             pos, t, dt) /
                     density;
}
}  // namespace MaterialPropertyLib
