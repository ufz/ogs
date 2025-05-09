/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TransportPorosity.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{

template <int DisplacementDim>

void TransportPorosityModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    SolidCompressibilityData const& solid_compressibility_data,
    BishopsData const& bishops_data,
    PrevState<BishopsData> const& bishops_data_prev,
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    PorosityData const& poro_data,
    ProcessLib::ConstitutiveRelations::MechanicalStrainData<
        DisplacementDim> const& eps_m_data,
    PrevState<ProcessLib::ConstitutiveRelations::MechanicalStrainData<
        DisplacementDim>> const& eps_m_prev_data,
    PrevState<TransportPorosityData> const& transport_poro_data_prev,
    TransportPorosityData& transport_poro_data) const
{
    namespace MPL = MaterialPropertyLib;

    auto const& medium = media_data.medium;

    if (!medium.hasProperty(MPL::PropertyType::transport_porosity))
    {
        transport_poro_data.phi = poro_data.phi;
        return;
    }

    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    static constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<kelvin_vector_size>;
    // Used in
    // MaterialLib/MPL/Properties/PermeabilityOrthotropicPowerLaw.cpp
    variables_prev.transport_porosity = transport_poro_data_prev->phi;

    // Used in
    // MaterialLib/MPL/Properties/TransportPorosityFromMassBalance.cpp
    variables.grain_compressibility = solid_compressibility_data.beta_SR;
    // Set volumetric mechanical strain rate for the general case without
    // swelling.
    variables.volumetric_mechanical_strain =
        Invariants::trace(eps_m_data.eps_m);
    variables_prev.volumetric_mechanical_strain =
        Invariants::trace(eps_m_prev_data->eps_m);
    variables.effective_pore_pressure =
        -bishops_data.chi_S_L * p_cap_data.p_cap;
    variables.porosity = poro_data.phi;

    // Used in MaterialLib/MPL/Properties/PorosityFromMassBalance.cpp
    // and MaterialLib/MPL/Properties/TransportPorosityFromMassBalance.cpp
    variables_prev.effective_pore_pressure =
        -bishops_data_prev->chi_S_L * p_cap_data.p_cap_prev;

    transport_poro_data.phi =
        medium.property(MPL::PropertyType::transport_porosity)
            .template value<double>(variables, variables_prev, x_t.x, x_t.t,
                                    x_t.dt);
}

template struct TransportPorosityModel<2>;
template struct TransportPorosityModel<3>;
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
