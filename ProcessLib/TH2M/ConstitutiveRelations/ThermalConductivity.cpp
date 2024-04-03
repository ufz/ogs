/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ThermalConductivity.h"

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void ThermalConductivityModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    TemperatureData const& T_data, PorosityData const& porosity_data,
    PorosityDerivativeData const& porosity_d_data,
    SaturationData const& S_L_data, SaturationDataDeriv const& dS_L_dp_cap,
    ThermalConductivityData<DisplacementDim>& thermal_conductivity_data) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables.temperature = T_data.T;
    variables.porosity = porosity_data.phi;
    variables.liquid_saturation = S_L_data.S_L;

    auto const& mpl_thermal_conductivity =
        media_data.medium[MPL::PropertyType::thermal_conductivity];

    thermal_conductivity_data.lambda = MPL::formEigenTensor<DisplacementDim>(
        mpl_thermal_conductivity.value(variables, x_t.x, x_t.t, x_t.dt));

    // Derivatives computed here and not in the MPL property because various
    // derivatives are not available in the VariableArray.

    auto const lambdaGR =
        media_data.gas.hasProperty(MPL::PropertyType::thermal_conductivity)
            ? MPL::formEigenTensor<DisplacementDim>(
                  media_data.gas[MPL::PropertyType::thermal_conductivity].value(
                      variables, x_t.x, x_t.t, x_t.dt))
            : MPL::formEigenTensor<DisplacementDim>(0.);

    auto const dlambda_GR_dT =
        media_data.gas.hasProperty(MPL::PropertyType::thermal_conductivity)
            ? MPL::formEigenTensor<DisplacementDim>(
                  media_data.gas[MPL::PropertyType::thermal_conductivity]
                      .dValue(variables, MPL::Variable::temperature, x_t.x,
                              x_t.t, x_t.dt))
            : MPL::formEigenTensor<DisplacementDim>(0.);

    auto const lambdaLR =
        media_data.liquid.hasProperty(MPL::PropertyType::thermal_conductivity)
            ? MPL::formEigenTensor<DisplacementDim>(
                  media_data.liquid[MPL::PropertyType::thermal_conductivity]
                      .value(variables, x_t.x, x_t.t, x_t.dt))
            : MPL::formEigenTensor<DisplacementDim>(0.);

    auto const dlambda_LR_dT =
        media_data.liquid.hasProperty(MPL::PropertyType::thermal_conductivity)
            ? MPL::formEigenTensor<DisplacementDim>(
                  media_data.liquid[MPL::PropertyType::thermal_conductivity]
                      .dValue(variables, MPL::Variable::temperature, x_t.x,
                              x_t.t, x_t.dt))
            : MPL::formEigenTensor<DisplacementDim>(0.);

    auto const lambdaSR =
        media_data.solid.hasProperty(MPL::PropertyType::thermal_conductivity)
            ? MPL::formEigenTensor<DisplacementDim>(
                  media_data.solid[MPL::PropertyType::thermal_conductivity]
                      .value(variables, x_t.x, x_t.t, x_t.dt))
            : MPL::formEigenTensor<DisplacementDim>(0.);

    auto const dlambda_SR_dT =
        media_data.solid.hasProperty(MPL::PropertyType::thermal_conductivity)
            ? MPL::formEigenTensor<DisplacementDim>(
                  media_data.solid[MPL::PropertyType::thermal_conductivity]
                      .dValue(variables, MPL::Variable::temperature, x_t.x,
                              x_t.t, x_t.dt))
            : MPL::formEigenTensor<DisplacementDim>(0.);

    // dphi_G_dp_GR = -ds_L_dp_GR * phi = 0;
    double const dphi_G_dp_cap = -dS_L_dp_cap() * porosity_data.phi;
    // dphi_L_dp_GR = ds_L_dp_GR * phi = 0;
    double const dphi_L_dp_cap = dS_L_dp_cap() * porosity_data.phi;

    thermal_conductivity_data.dlambda_dp_cap =
        dphi_G_dp_cap * lambdaGR + dphi_L_dp_cap * lambdaLR;

    double const phi_L = S_L_data.S_L * porosity_data.phi;
    double const phi_G = (1. - S_L_data.S_L) * porosity_data.phi;
    double const phi_S = 1. - porosity_data.phi;

    thermal_conductivity_data.dlambda_dT =
        phi_G * dlambda_GR_dT + phi_L * dlambda_LR_dT + phi_S * dlambda_SR_dT -
        porosity_d_data.dphi_dT * lambdaSR;
}
template struct ThermalConductivityModel<2>;
template struct ThermalConductivityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
