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
    SaturationData const& S_L_data,
    TemperatureData<DisplacementDim> const& T_data,
    PorosityData const& poro_data, LiquidViscosityData const& mu_L_data,
    PorosityData& transport_poro_data,
    PorosityData const& transport_poro_data_prev,
    SolidMechanicsDataStateless<DisplacementDim> const& s_mech_data,
    PermeabilityData<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;

    auto const& medium = media_data.medium;

    MPL::VariableArray variables;
    variables.liquid_saturation = S_L_data.S_L;
    variables.temperature = T_data.T;
    MPL::VariableArray variables_prev;

    if (medium.hasProperty(MPL::PropertyType::transport_porosity))
    {
        // Used in
        // MaterialLib/MPL/Properties/PermeabilityOrthotropicPowerLaw.cpp
        variables_prev.transport_porosity = transport_poro_data_prev.phi;

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
