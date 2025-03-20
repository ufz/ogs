/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SolidDensity.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void SolidDensityModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t,
    MediaData const& media_data,
    PorosityData const& poro_data,
    TemperatureData<DisplacementDim> const& T_data,
    EffectiveStressData<DisplacementDim>& sigma_eff_data,
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    BishopsData const& bishops_data,
    SolidDensityData& out) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables.temperature = T_data.T;

    double const p_FR = -bishops_data.chi_S_L * p_cap_data.p_cap;

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;
    variables.solid_grain_pressure =
        p_FR -
        sigma_eff_data.sigma_eff.dot(identity2) / (3 * (1 - poro_data.phi));

    out.rho_SR = media_data.solid.property(MPL::PropertyType::density)
                     .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
    out.dry_density_solid =
        (1 - poro_data.phi) * out.rho_SR;  // TODO only for output
}

template struct SolidDensityModel<2>;
template struct SolidDensityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
