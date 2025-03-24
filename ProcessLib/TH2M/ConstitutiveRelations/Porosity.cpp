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
    SpaceTimeData const& x_t, MediaData const& media_data,
    SaturationData const& S_L_data,
    PrevState<SaturationData> const& S_L_prev_data,
    CapillaryPressureData const& p_cap, GasPressureData const& p_GR,
    BishopsData const& chi_S_L, PrevState<BishopsData> const& chi_S_L_prev,
    BiotData const& biot, SolidCompressibilityData const& solid_compressibility,
    StrainData<DisplacementDim> const& eps_data,
    PrevState<StrainData<DisplacementDim>> const& eps_prev_data,
    PrevState<PorosityData> const& porosity_prev_data,
    PorosityData& porosity_data) const
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    MaterialPropertyLib::VariableArray variables;
    MaterialPropertyLib::VariableArray variables_prev;

    variables.grain_compressibility = solid_compressibility();

    variables.liquid_saturation = S_L_data.S_L;
    variables_prev.liquid_saturation = S_L_prev_data->S_L;

    /*auto const& mpl_porosity =
        media_data.medium[MaterialPropertyLib::PropertyType::porosity];

    double const phi_0 =
        mpl_porosity.template value<double>(variables, x_t.x, x_t.t, x_t.dt);
    */

    variables.effective_pore_pressure =
        (1 - chi_S_L.chi_S_L) * p_GR.pG + chi_S_L.chi_S_L * (p_GR.pG - p_cap.pCap);

    // Used in MaterialLib/MPL/Properties/PorosityFromMassBalance.cpp
    // and MaterialLib/MPL/Properties/TransportPorosityFromMassBalance.cpp
    variables_prev.effective_pore_pressure =
        chi_S_L_prev->chi_S_L *
        ((1 - S_L_prev_data->S_L) * p_GR.pG_prev +
         S_L_prev_data->S_L * (p_GR.pG_prev - p_cap.pCap_prev));

    variables.volumetric_strain = Invariants::trace(eps_data.eps);
    variables_prev.volumetric_strain = Invariants::trace(eps_prev_data->eps);

    /*double const phi_prev = porosity_prev_data->phi;
    if (std::isnan(phi_prev) == true)
    {
        phi_prev = phi_0;
    }*/

    // double const phi_S =
    //     (1. - phi_0) *
    //    (1. - biot() * div_u);
    /*double const phi_S =
        (1. - phi_0) *
        (1. + s_therm_exp_data.thermal_volume_strain - biot() * div_u);

    porosity_data.phi = 1. - phi_S;*/
    variables_prev.porosity = porosity_prev_data->phi;
    porosity_data.phi =
        media_data.medium.property(MaterialPropertyLib::PropertyType::porosity)
            .template value<double>(variables, variables_prev, x_t.x, x_t.t,
                                    x_t.dt);
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

    /*auto const& mpl_porosity =
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

    porosity_d_data.dphi_L_dp_cap = dS_L_dp_cap() * porosity_data.phi;*/
    porosity_d_data.dphi_dT = 0.0;
    porosity_d_data.dphi_L_dp_cap = 0.0;
}

template struct PorosityModelNonConstantSolidPhaseVolumeFraction<2>;
template struct PorosityModelNonConstantSolidPhaseVolumeFraction<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
