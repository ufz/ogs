// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SolidThermalExpansion.h"

#include "MaterialLib/MPL/Utils/FormKelvinVector.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void SolidThermalExpansionModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    SolidThermalExpansionData<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;

    // Consider also anisotropic thermal expansion.
    out.solid_linear_thermal_expansivity_vector =
        MPL::formKelvinVector<DisplacementDim>(
            media_data.solid.property(MPL::PropertyType::thermal_expansivity)
                .value(variables, x_t.x, x_t.t, x_t.dt));
}

template struct SolidThermalExpansionModel<2>;
template struct SolidThermalExpansionModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
