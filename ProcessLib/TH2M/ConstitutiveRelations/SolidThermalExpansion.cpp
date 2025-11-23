/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SolidThermalExpansion.h"

#include "MaterialLib/MPL/Utils/FormKelvinVector.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void SolidThermalExpansionModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    TemperatureData const& T_data, ReferenceTemperatureData T0,
    SolidThermalExpansionData<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;

    // Consider also anisotropic thermal expansion.
    out.solid_linear_thermal_expansivity_vector =
        MPL::formKelvinVector<DisplacementDim>(
            media_data.thermal_expansivity_solid.value(variables, x_t.x, x_t.t,
                                                       x_t.dt));

    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    out.beta_T_SR =
        Invariants::trace(out.solid_linear_thermal_expansivity_vector);

    double const delta_T = T_data.T - T0();
    out.thermal_volume_strain = out.beta_T_SR * delta_T;
}

template struct SolidThermalExpansionModel<2>;
template struct SolidThermalExpansionModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
