/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PermeabilityModel.h"

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void PermeabilityModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    SaturationData const& S_L_data, GasPressureData const& p_g,
    CapillaryPressureData const& p_cap, TemperatureData const& T_data,
    TransportPorosityData const& transport_poro_data,
    TotalStressData<DisplacementDim> const& total_stress_data,
    MechanicalStrainData<DisplacementDim> const& mechanical_strain_data,
    StrainData<DisplacementDim> const& eps_data,
    EquivalentPlasticStrainData const& equivalent_plastic_strain,
    PermeabilityData<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;

    static constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<kelvin_vector_size>;

    auto const& medium = media_data.medium;

    MPL::VariableArray variables;
    variables.liquid_saturation = S_L_data.S_L;
    variables.temperature = T_data.T;
    variables.capillary_pressure = p_cap.pCap;
    variables.gas_phase_pressure = p_g.pG;
    variables.transport_porosity = transport_poro_data.phi;

    out.k_rel_G =
        medium
            .property(MPL::PropertyType::relative_permeability_nonwetting_phase)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    out.dk_rel_G_dS_L =
        medium[MPL::PropertyType::relative_permeability_nonwetting_phase]
            .template dValue<double>(variables,
                                     MPL::Variable::liquid_saturation, x_t.x,
                                     x_t.t, x_t.dt);

    out.k_rel_L = medium.property(MPL::PropertyType::relative_permeability)
                      .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    out.dk_rel_L_dS_L = medium[MPL::PropertyType::relative_permeability]
                            .template dValue<double>(
                                variables, MPL::Variable::liquid_saturation,
                                x_t.x, x_t.t, x_t.dt);

    // For stress dependent permeability.
    variables.total_stress =
        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
            total_stress_data.sigma_total);

    variables.mechanical_strain.emplace<KelvinVector<DisplacementDim>>(
        mechanical_strain_data.eps_m);

    variables.equivalent_plastic_strain = equivalent_plastic_strain();

    variables.volumetric_strain = Invariants::trace(eps_data.eps);

    out.Ki = MPL::formEigenTensor<DisplacementDim>(
        medium.property(MPL::PropertyType::permeability)
            .value(variables, x_t.x, x_t.t, x_t.dt));
}

template struct PermeabilityModel<2>;
template struct PermeabilityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
