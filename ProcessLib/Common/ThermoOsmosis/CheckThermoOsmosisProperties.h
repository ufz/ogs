// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <memory>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Phase.h"
#include "MaterialLib/MPL/PropertyType.h"

namespace ProcessLib
{
/// Checks the thermo-osmosis parametrisation of the given medium: the
/// thermo-osmosis properties must be defined on the medium and not on the
/// solid phase, and the two mutually exclusive parametrisations
/// \c thermal_osmosis_coefficient and \c thermal_osmosis_permeability must not
/// be defined at the same time.
///
/// The properties themselves are optional; a medium without any of them has no
/// thermo-osmosis. See getThermoOsmoticCoefficient() for their meaning.
inline void checkThermoOsmosisProperties(
    MaterialPropertyLib::Medium const& medium)
{
    auto const solid_phase =
        getOptionalPhase(medium, MaterialPropertyLib::PhaseName::Solid);
    if (solid_phase)
    {
        // thermal_osmosis_coefficient used to live on the solid phase;
        // thermal_osmosis_permeability never did, but a project file
        // hand-migrated from the former can put it there just as easily.
        for (auto const property_type :
             {MaterialPropertyLib::PropertyType::thermal_osmosis_coefficient,
              MaterialPropertyLib::PropertyType::thermal_osmosis_permeability})
        {
            if (solid_phase->hasProperty(property_type))
            {
                OGS_FATAL(
                    "{:s} is defined on the solid phase of {:s}, but is read "
                    "from the medium. Move the property from the solid phase "
                    "to the medium, or use "
                    "scripts/dev/move_thermal_osmosis.py to migrate the "
                    "project file.",
                    MaterialPropertyLib::property_enum_to_string[property_type],
                    medium.description());
            }
        }
    }

    if (medium.hasProperty(
            MaterialPropertyLib::PropertyType::thermal_osmosis_permeability) &&
        medium.hasProperty(
            MaterialPropertyLib::PropertyType::thermal_osmosis_coefficient))
    {
        OGS_FATAL(
            "Thermo-osmosis permeability and coefficient cannot be defined at "
            "the same time in {:s}.",
            medium.description());
    }
}

/// Checks the thermo-osmosis parametrisation of each of the given media, see
/// checkThermoOsmosisProperties(MaterialPropertyLib::Medium const&).
inline void checkThermoOsmosisProperties(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    for (auto const& m : media)
    {
        checkThermoOsmosisProperties(*m.second);
    }
}
}  // namespace ProcessLib
