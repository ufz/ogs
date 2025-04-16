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

template <int DisplacementDim>
void PorosityModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    SaturationData const& S_L_data,
    PrevState<SaturationData> const& S_L_prev_data,
    CapillaryPressureData const& p_cap, GasPressureData const& p_GR,
    BishopsData const& chi_S_L, PrevState<BishopsData> const& chi_S_L_prev,
    SolidCompressibilityData const& solid_compressibility,
    StrainData<DisplacementDim> const& eps_data,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const&& eps_prev,
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

    variables.effective_pore_pressure =
        (1 - chi_S_L.chi_S_L) * p_GR.pG +
        chi_S_L.chi_S_L * (p_GR.pG - p_cap.pCap);

    // Used in MaterialLib/MPL/Properties/PorosityFromMassBalance.cpp
    // and MaterialLib/MPL/Properties/TransportPorosityFromMassBalance.cpp
    variables_prev.effective_pore_pressure =
        (1 - chi_S_L_prev->chi_S_L) * p_GR.pG_prev +
        chi_S_L_prev->chi_S_L * (p_GR.pG_prev - p_cap.pCap_prev);

    variables.volumetric_strain = Invariants::trace(eps_data.eps);
    variables_prev.volumetric_strain = Invariants::trace(eps_prev);

    variables_prev.porosity = porosity_prev_data->phi;
    porosity_data.phi =
        media_data.medium.property(MaterialPropertyLib::PropertyType::porosity)
            .template value<double>(variables, variables_prev, x_t.x, x_t.t,
                                    x_t.dt);
}

template <int DisplacementDim>
void PorosityModel<DisplacementDim>::dEval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    SaturationData const& S_L_data,
    PrevState<SaturationData> const& S_L_prev_data,
    CapillaryPressureData const& p_cap, GasPressureData const& p_GR,
    BishopsData const& chi_S_L, PrevState<BishopsData> const& chi_S_L_prev,
    SolidCompressibilityData const& solid_compressibility,
    StrainData<DisplacementDim> const& eps_data,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const&& eps_prev,
    PrevState<PorosityData> const& porosity_prev_data,
    PorosityDerivativeData& porosity_d_data) const
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

    variables.volumetric_strain = Invariants::trace(eps_data.eps);
    variables_prev.volumetric_strain = Invariants::trace(eps_prev);

    variables_prev.porosity = porosity_prev_data->phi;

    auto const& mpl_porosity =
        media_data.medium[MaterialPropertyLib::PropertyType::porosity];

    porosity_d_data.dphi_dT = mpl_porosity.template dValue<double>(
        variables, variables_prev, MaterialPropertyLib::Variable::temperature,
        x_t.x, x_t.t, x_t.dt);

    porosity_d_data.dphi_L_dp_cap = mpl_porosity.template dValue<double>(
        variables, variables_prev,
        MaterialPropertyLib::Variable::capillary_pressure, x_t.x, x_t.t,
        x_t.dt);
}

template struct PorosityModel<2>;
template struct PorosityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
