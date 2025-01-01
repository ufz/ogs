/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
    SaturationData const& S_L_data,
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
}

template <int DisplacementDim>
void ThermalConductivityModel<DisplacementDim>::dEval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    TemperatureData const& T_data, PorosityData const& porosity_data,
    PorosityDerivativeData const& porosity_d_data,
    SaturationData const& S_L_data,
    ThermalConductivityDerivativeData<DisplacementDim>&
        thermal_conductivity_d_data) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables.temperature = T_data.T;
    variables.porosity = porosity_data.phi;
    variables.liquid_saturation = S_L_data.S_L;

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

    thermal_conductivity_d_data.dlambda_dp_cap =
        -porosity_d_data.dphi_L_dp_cap * lambdaGR +
        porosity_d_data.dphi_L_dp_cap * lambdaLR;

    double const phi_L = S_L_data.S_L * porosity_data.phi;
    double const phi_G = (1. - S_L_data.S_L) * porosity_data.phi;
    double const phi_S = 1. - porosity_data.phi;

    // Assuming dS_L/dT = 0, then:
    // dphi_G_dT = -dS_L/dT * phi + (1 - S_L) * dphi_dT = (1 - S_L) * dphi_dT
    // dphi_L_dT =  dS_L/dT * phi +      S_L  * dphi_dT        S_L  * dphi_dT
    // dphi_S_dT =                             -dphi_dT              -dphi_dT
    thermal_conductivity_d_data.dlambda_dT =
        (1 - S_L_data.S_L) * porosity_d_data.dphi_dT * lambdaGR +
        phi_G * dlambda_GR_dT +
        S_L_data.S_L * porosity_d_data.dphi_dT * lambdaLR +
        +phi_L * dlambda_LR_dT - porosity_d_data.dphi_dT * lambdaSR +
        phi_S * dlambda_SR_dT;
}

template struct ThermalConductivityModel<2>;
template struct ThermalConductivityModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
