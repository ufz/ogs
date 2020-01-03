/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on August 16, 2019, 3:40 PM
 */

#include "GetThermalExpansivity.h"

#include "MaterialLib/MPL/Phase.h"

namespace MaterialPropertyLib
{
class Phase;

double getThermalExpansivity(Phase const& phase, VariableArray const& vars,
                             const double density,
                             ParameterLib::SpatialPosition const& pos,
                             double const t)
{
    auto const thermal_expansivity_ptr =
        &phase.property(MaterialPropertyLib::PropertyType::thermal_expansivity);

    // The thermal expansivity is explicitly given in the project file.
    if (thermal_expansivity_ptr)
    {
        return (*thermal_expansivity_ptr).template value<double>(vars, pos, t);
    }

    // The thermal expansivity calculated by the density model directly.
    return (density == 0.0)
               ? 0.0
               : -phase.property(MaterialPropertyLib::PropertyType::density)
                         .template dValue<double>(
                             vars, MaterialPropertyLib::Variable::temperature,
                             pos, t) /
                     density;
}
}  // namespace MaterialPropertyLib
