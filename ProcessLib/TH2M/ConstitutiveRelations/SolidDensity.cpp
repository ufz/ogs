/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SolidDensity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{

template <int DisplacementDim>
void SolidDensityModel<DisplacementDim>::eval(SpaceTimeData const& x_t,
                             MediaData const& media_data,
                             TemperatureData const& T_data,
                             ProcessLib::ConstitutiveRelations::EffectiveStressData<
                                    DisplacementDim> const& sigma_eff_data,
                             CapillaryPressureData const& p_cap,
                             GasPressureData const& p_GR,
                             BishopsData const& chi_S_L,
                             PorosityData const& poro_data,
                             SolidDensityData& solid_density_data) const
{
    MaterialPropertyLib::VariableArray variables;
    variables.temperature = T_data.T;

    double const p_FR = (1 - chi_S_L.chi_S_L) * p_GR.pG + chi_S_L.chi_S_L * (p_GR.pG - p_cap.pCap);
    
    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    variables.solid_grain_pressure =
            p_FR -
            sigma_eff_data.sigma_eff.dot(identity2) / (3 * (1 - poro_data.phi));
    
    auto const& mpl_solid_density =
        media_data.solid[MaterialPropertyLib::PropertyType::density];

    solid_density_data.rho_SR = mpl_solid_density.template value<double>(
        variables, x_t.x, x_t.t, x_t.dt);
}

template <int DisplacementDim>
void SolidDensityModel<DisplacementDim>::dEval(
    SpaceTimeData const& x_t,
    MediaData const& media_data,
    TemperatureData const& T_data,
    ProcessLib::ConstitutiveRelations::EffectiveStressData<
                                    DisplacementDim> const& sigma_eff_data,
    CapillaryPressureData const& p_cap,
    GasPressureData const& p_GR,
    BishopsData const& chi_S_L,
    PorosityData const& poro_data,
    SolidDensityDerivativeData& solid_density_d_data) const
{
    MaterialPropertyLib::VariableArray variables;
    variables.temperature = T_data.T;

    double const p_FR = (1 - chi_S_L.chi_S_L) * p_GR.pG + chi_S_L.chi_S_L * (p_GR.pG - p_cap.pCap);
    
    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    variables.solid_grain_pressure =
            p_FR -
            sigma_eff_data.sigma_eff.dot(identity2) / (3 * (1 - poro_data.phi));
    

    auto const& mpl_solid_density =
        media_data.solid[MaterialPropertyLib::PropertyType::density];

    solid_density_d_data.drho_SR_dT = mpl_solid_density.template dValue<double>(
        variables, MaterialPropertyLib::Variable::temperature, x_t.x, x_t.t,
        x_t.dt);
}

template struct SolidDensityModel<2>;
template struct SolidDensityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
