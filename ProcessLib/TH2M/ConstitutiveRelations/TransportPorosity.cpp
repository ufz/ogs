/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TransportPorosity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{

template <int DisplacementDim>
void TransportPorosityModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    SaturationData const& S_L_data,
    PrevState<SaturationData> const& S_L_prev_data,
    CapillaryPressureData const& p_cap, GasPressureData const& p_GR,
    BishopsData const& chi_S_L, PrevState<BishopsData> const& chi_S_L_prev,
    SolidCompressibilityData const& solid_compressibility,
    MechanicalStrainData<DisplacementDim> const& eps_m_data,
    PrevState<MechanicalStrainData<DisplacementDim>> const& eps_m_prev_data,
    PrevState<TransportPorosityData> const& transport_porosity_prev_data,
    PorosityData const& poro_data,
    TransportPorosityData& transport_porosity_data) const
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    MaterialPropertyLib::VariableArray variables;
    MaterialPropertyLib::VariableArray variables_prev;

    variables.grain_compressibility = solid_compressibility();

    variables.liquid_saturation = S_L_data.S_L;
    variables_prev.liquid_saturation = S_L_prev_data->S_L;

    variables.effective_pore_pressure =
        (1 - chi_S_L.chi_S_L) * p_GR.pG +
        chi_S_L.chi_S_L * (p_GR.pG - p_cap.pCap);

    // Used in MaterialLib/MPL/Properties/PorosityFromMassBalance.cpp
    // and MaterialLib/MPL/Properties/TransportPorosityFromMassBalance.cpp
    variables_prev.effective_pore_pressure =
        (1 - chi_S_L_prev->chi_S_L) * p_GR.pG_prev +
        chi_S_L_prev->chi_S_L * (p_GR.pG_prev - p_cap.pCap_prev);

    variables.volumetric_strain = Invariants::trace(eps_m_data.eps_m);
    variables_prev.volumetric_strain =
        Invariants::trace(eps_m_prev_data->eps_m);

    variables_prev.transport_porosity = transport_porosity_prev_data->phi;
    variables.porosity = poro_data.phi;
    transport_porosity_data.phi =
        media_data.medium
            .property(MaterialPropertyLib::PropertyType::transport_porosity)
            .template value<double>(variables, variables_prev, x_t.x, x_t.t,
                                    x_t.dt);
}

template <int DisplacementDim>
void TransportPorosityModel<DisplacementDim>::dEval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    SaturationData const& S_L_data,
    PrevState<SaturationData> const& S_L_prev_data,
    CapillaryPressureData const& p_cap, GasPressureData const& p_GR,
    BishopsData const& chi_S_L, PrevState<BishopsData> const& chi_S_L_prev,
    SolidCompressibilityData const& solid_compressibility,
    MechanicalStrainData<DisplacementDim> const& eps_m_data,
    PrevState<MechanicalStrainData<DisplacementDim>> const& eps_m_prev_data,
    PrevState<TransportPorosityData> const& transport_porosity_prev_data,
    PorosityData const& poro_data,
    TransportPorosityDerivativeData& transport_porosity_d_data) const
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    MaterialPropertyLib::VariableArray variables;
    MaterialPropertyLib::VariableArray variables_prev;

    variables.grain_compressibility = solid_compressibility();

    variables.liquid_saturation = S_L_data.S_L;
    variables_prev.liquid_saturation = S_L_prev_data->S_L;

    variables.effective_pore_pressure =
        (1 - chi_S_L.chi_S_L) * p_GR.pG +
        chi_S_L.chi_S_L * (p_GR.pG - p_cap.pCap);

    // Used in MaterialLib/MPL/Properties/PorosityFromMassBalance.cpp
    // and MaterialLib/MPL/Properties/TransportPorosityFromMassBalance.cpp
    variables_prev.effective_pore_pressure =
        (1 - chi_S_L_prev->chi_S_L) * p_GR.pG_prev +
        chi_S_L_prev->chi_S_L * (p_GR.pG_prev - p_cap.pCap_prev);

    variables.volumetric_strain = Invariants::trace(eps_m_data.eps_m);
    variables_prev.volumetric_strain =
        Invariants::trace(eps_m_prev_data->eps_m);

    variables_prev.transport_porosity = transport_porosity_prev_data->phi;
    variables.porosity = poro_data.phi;

    auto const& mpl_transport_porosity =
        media_data
            .medium[MaterialPropertyLib::PropertyType::transport_porosity];

    transport_porosity_d_data.dphi_dT =
        mpl_transport_porosity.template dValue<double>(
            variables, MaterialPropertyLib::Variable::temperature, x_t.x, x_t.t,
            x_t.dt);

    transport_porosity_d_data.dphi_L_dp_cap =
        mpl_transport_porosity.template dValue<double>(
            variables, MaterialPropertyLib::Variable::capillary_pressure, x_t.x,
            x_t.t, x_t.dt);
}

template struct TransportPorosityModel<2>;
template struct TransportPorosityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
