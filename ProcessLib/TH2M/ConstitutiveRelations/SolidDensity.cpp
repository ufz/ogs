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

void SolidDensityModel::eval(SpaceTimeData const& x_t,
                             MediaData const& media_data,
                             TemperatureData const& T_data,
                             SolidDensityData& solid_density_data) const
{
    MaterialPropertyLib::VariableArray variables;
    variables.temperature = T_data.T;

    auto const& mpl_solid_density =
        media_data.solid[MaterialPropertyLib::PropertyType::density];

    solid_density_data.rho_SR = mpl_solid_density.template value<double>(
        variables, x_t.x, x_t.t, x_t.dt);
}

void SolidDensityModel::dEval(
    SpaceTimeData const& x_t,
    MediaData const& media_data,
    TemperatureData const& T_data,
    SolidDensityDerivativeData& solid_density_d_data) const
{
    MaterialPropertyLib::VariableArray variables;
    variables.temperature = T_data.T;

    auto const& mpl_solid_density =
        media_data.solid[MaterialPropertyLib::PropertyType::density];

    solid_density_d_data.drho_SR_dT = mpl_solid_density.template dValue<double>(
        variables, MaterialPropertyLib::Variable::temperature, x_t.x, x_t.t,
        x_t.dt);
}

template <int DisplacementDim>
void SolidDensityModelNonConstantSolidPhaseVolumeFraction<DisplacementDim>::
    eval(SpaceTimeData const& x_t,
         MediaData const& media_data,
         TemperatureData const& T_data,
         BiotData const& biot,
         StrainData<DisplacementDim> const& strain_data,
         SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
         SolidDensityData& solid_density_data) const
{
    MaterialPropertyLib::VariableArray variables;
    variables.temperature = T_data.T;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    double const div_u = Invariants::trace(strain_data.eps);

    auto const& mpl_solid_density =
        media_data.solid[MaterialPropertyLib::PropertyType::density];

    auto const rho_ref_SR = mpl_solid_density.template value<double>(
        variables, x_t.x, x_t.t, x_t.dt);

    solid_density_data.rho_SR =
        rho_ref_SR *
        (1. - s_therm_exp_data.thermal_volume_strain + (biot() - 1.) * div_u);
}

template <int DisplacementDim>
void SolidDensityModelNonConstantSolidPhaseVolumeFraction<DisplacementDim>::
    dEval(SpaceTimeData const& x_t,
          MediaData const& media_data,
          TemperatureData const& T_data,
          BiotData const& biot,
          StrainData<DisplacementDim> const& strain_data,
          SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
          SolidDensityDerivativeData& solid_density_d_data) const
{
    MaterialPropertyLib::VariableArray variables;
    variables.temperature = T_data.T;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    double const div_u = Invariants::trace(strain_data.eps);

    auto const& mpl_solid_density =
        media_data.solid[MaterialPropertyLib::PropertyType::density];

    auto const rho_ref_SR = mpl_solid_density.template value<double>(
        variables, x_t.x, x_t.t, x_t.dt);

    solid_density_d_data.drho_SR_dT =
        mpl_solid_density.template dValue<double>(
            variables, MaterialPropertyLib::Variable::temperature, x_t.x, x_t.t,
            x_t.dt) *
            (1. - s_therm_exp_data.thermal_volume_strain +
             (biot() - 1.) * div_u) -
        rho_ref_SR * s_therm_exp_data.beta_T_SR;
}

template struct SolidDensityModelNonConstantSolidPhaseVolumeFraction<2>;
template struct SolidDensityModelNonConstantSolidPhaseVolumeFraction<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
