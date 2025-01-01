/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Porosity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{

void PorosityModel::eval(SpaceTimeData const& x_t,
                         MediaData const& media_data,
                         PorosityData& porosity_data) const
{
    MaterialPropertyLib::VariableArray variables;

    auto const& mpl_porosity =
        media_data.medium[MaterialPropertyLib::PropertyType::porosity];

    porosity_data.phi =
        mpl_porosity.template value<double>(variables, x_t.x, x_t.t, x_t.dt);
}

void PorosityModel::dEval(SpaceTimeData const& x_t, MediaData const& media_data,
                          PorosityData const& porosity_data,
                          SaturationDataDeriv const& dS_L_dp_cap,
                          PorosityDerivativeData& porosity_d_data) const
{
    MaterialPropertyLib::VariableArray variables;

    auto const& mpl_porosity =
        media_data.medium[MaterialPropertyLib::PropertyType::porosity];

    porosity_d_data.dphi_dT = mpl_porosity.template dValue<double>(
        variables, MaterialPropertyLib::Variable::temperature, x_t.x, x_t.t,
        x_t.dt);

    porosity_d_data.dphi_L_dp_cap = dS_L_dp_cap() * porosity_data.phi;
}

template <int DisplacementDim>
void PorosityModelNonConstantSolidPhaseVolumeFraction<DisplacementDim>::eval(
    SpaceTimeData const& x_t,
    MediaData const& media_data,
    BiotData const& biot,
    StrainData<DisplacementDim> const& strain_data,
    SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
    PorosityData& porosity_data) const
{
    MaterialPropertyLib::VariableArray variables;

    auto const& mpl_porosity =
        media_data.medium[MaterialPropertyLib::PropertyType::porosity];

    double const phi_0 =
        mpl_porosity.template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    double const div_u = Invariants::trace(strain_data.eps);

    double const phi_S =
        (1. - phi_0) *
        (1. + s_therm_exp_data.thermal_volume_strain - biot() * div_u);

    porosity_data.phi = 1. - phi_S;
}

template <int DisplacementDim>
void PorosityModelNonConstantSolidPhaseVolumeFraction<DisplacementDim>::dEval(
    SpaceTimeData const& x_t,
    MediaData const& media_data,
    PorosityData const& porosity_data,
    SaturationDataDeriv const& dS_L_dp_cap,
    BiotData const& biot,
    SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
    StrainData<DisplacementDim> const& strain_data,
    PorosityDerivativeData& porosity_d_data) const
{
    MaterialPropertyLib::VariableArray variables;

    auto const& mpl_porosity =
        media_data.medium[MaterialPropertyLib::PropertyType::porosity];

    double const phi_0 =
        mpl_porosity.template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    double const div_u = Invariants::trace(strain_data.eps);

    auto const dphi_0_dT = mpl_porosity.template dValue<double>(
        variables, MaterialPropertyLib::Variable::temperature, x_t.x, x_t.t,
        x_t.dt);

    porosity_d_data.dphi_dT =
        dphi_0_dT *
            (1. + s_therm_exp_data.thermal_volume_strain - biot() * div_u) -
        (1. - phi_0) * s_therm_exp_data.beta_T_SR;

    porosity_d_data.dphi_L_dp_cap = dS_L_dp_cap() * porosity_data.phi;
}

template struct PorosityModelNonConstantSolidPhaseVolumeFraction<2>;
template struct PorosityModelNonConstantSolidPhaseVolumeFraction<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
