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
/// Checks that the liquid phase of the given medium has no
/// \c thermal_expansivity property.
///
/// The liquid thermal expansivity is computed from the temperature derivative
/// of the liquid density, \f$-\frac{1}{\rho_L}\frac{\partial \rho_L}{\partial
/// T}\f$. A \c thermal_expansivity property on the liquid phase is not read,
/// and is rejected to avoid silently ignoring the input.
inline void checkLiquidThermalExpansivity(
    MaterialPropertyLib::Medium const& medium)
{
    auto const liquid_phase =
        getOptionalPhase(medium, MaterialPropertyLib::PhaseName::AqueousLiquid);
    if (liquid_phase &&
        liquid_phase->hasProperty(
            MaterialPropertyLib::PropertyType::thermal_expansivity))
    {
        OGS_FATAL(
            "thermal_expansivity is defined on the AqueousLiquid phase of "
            "{:s}, but it is not used. The liquid thermal expansivity is "
            "computed from the temperature derivative of the liquid density. "
            "Remove the property and define a temperature dependent liquid "
            "density instead.",
            medium.description());
    }
}

/// Checks the liquid phase of each of the given media, see
/// checkLiquidThermalExpansivity(MaterialPropertyLib::Medium const&).
inline void checkLiquidThermalExpansivity(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    for (auto const& m : media)
    {
        checkLiquidThermalExpansivity(*m.second);
    }
}
}  // namespace ProcessLib
