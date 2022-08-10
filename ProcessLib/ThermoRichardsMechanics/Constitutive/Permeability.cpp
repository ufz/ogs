/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Permeability.h"

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void PermeabilityModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    SolidCompressibilityData const& solid_compressibility_data,
    SaturationData const& S_L_data, BishopsData const& bishops_data,
    BishopsData const& bishops_data_prev,
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    TemperatureData<DisplacementDim> const& T_data,
    PorosityData const& poro_data, LiquidViscosityData const& mu_L_data,
    TransportPorosityData& transport_poro_data,
    TransportPorosityData const& transport_poro_data_prev,
    SolidMechanicsDataStateless<DisplacementDim> const& s_mech_data,
    StrainData<DisplacementDim> const& eps_data,
    StrainData<DisplacementDim> const& eps_prev_data,
    PermeabilityData<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;

    auto const& medium = media_data.medium;

    MPL::VariableArray variables;
    variables.liquid_saturation = S_L_data.S_L;
    variables.temperature = T_data.T;
    variables.capillary_pressure = p_cap_data.p_cap;
    MPL::VariableArray variables_prev;

    if (medium.hasProperty(MPL::PropertyType::transport_porosity))
    {
        static constexpr int kelvin_vector_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
        using Invariants =
            MathLib::KelvinVector::Invariants<kelvin_vector_size>;
        // Used in
        // MaterialLib/MPL/Properties/PermeabilityOrthotropicPowerLaw.cpp
        variables_prev.transport_porosity = transport_poro_data_prev.phi;

        // Used in
        // MaterialLib/MPL/Properties/TransportPorosityFromMassBalance.cpp
        variables.grain_compressibility = solid_compressibility_data.beta_SR;
        // Set volumetric strain rate for the general case without swelling.
        variables.volumetric_strain = Invariants::trace(eps_data.eps);
        variables_prev.volumetric_strain = Invariants::trace(eps_prev_data.eps);
        variables.effective_pore_pressure =
            -bishops_data.chi_S_L * p_cap_data.p_cap;
        variables.porosity = poro_data.phi;

        // Used in MaterialLib/MPL/Properties/PorosityFromMassBalance.cpp
        // and MaterialLib/MPL/Properties/TransportPorosityFromMassBalance.cpp
        variables_prev.effective_pore_pressure =
            -bishops_data_prev.chi_S_L *
            (p_cap_data.p_cap - p_cap_data.p_cap_dot * x_t.dt);

        transport_poro_data.phi =
            medium.property(MPL::PropertyType::transport_porosity)
                .template value<double>(variables, variables_prev, x_t.x, x_t.t,
                                        x_t.dt);
        variables.transport_porosity = transport_poro_data.phi;
    }
    else
    {
        variables.transport_porosity = poro_data.phi;
    }

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
            s_mech_data.sigma_total));

    variables.equivalent_plastic_strain = s_mech_data.equivalent_plastic_strain;

    auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
        medium.property(MPL::PropertyType::permeability)
            .value(variables, x_t.x, x_t.t, x_t.dt));

    out.Ki_over_mu = K_intrinsic / mu_L_data.viscosity;
}

template struct PermeabilityModel<2>;
template struct PermeabilityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
