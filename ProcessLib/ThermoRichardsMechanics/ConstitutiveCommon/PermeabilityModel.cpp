/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PermeabilityModel.h"

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void PermeabilityModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    SaturationData const& S_L_data,
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    TemperatureData<DisplacementDim> const& T_data,
    TransportPorosityData const& transport_poro_data,
    TotalStressData<DisplacementDim> const& total_stress_data,
    StrainData<DisplacementDim> const& eps_data,
    EquivalentPlasticStrainData const& equiv_plast_strain_data,
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
    variables.capillary_pressure = p_cap_data.p_cap;
    variables.transport_porosity = transport_poro_data.phi;

    out.k_rel = medium.property(MPL::PropertyType::relative_permeability)
                    .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
    out.dk_rel_dS_L =
        medium.property(MPL::PropertyType::relative_permeability)
            .template dValue<double>(variables,
                                     MPL::Variable::liquid_saturation,
                                     x_t.x,
                                     x_t.t,
                                     x_t.dt);

    // For stress dependent permeability.
    using SymmetricTensor =
        KelvinVector<DisplacementDim>;  // same data type, but different
                                        // semantics
    variables.total_stress.emplace<SymmetricTensor>(
        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
            total_stress_data.sigma_total));

    variables.equivalent_plastic_strain =
        equiv_plast_strain_data.equivalent_plastic_strain;

    variables.volumetric_strain = Invariants::trace(eps_data.eps);

    auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
        medium.property(MPL::PropertyType::permeability)
            .value(variables, x_t.x, x_t.t, x_t.dt));

    out.Ki = K_intrinsic;
}

template struct PermeabilityModel<2>;
template struct PermeabilityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
