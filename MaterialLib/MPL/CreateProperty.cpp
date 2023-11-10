/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateProperty.h"

#include <boost/algorithm/string/predicate.hpp>
#include <string>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "Component.h"
#include "Medium.h"
#include "Phase.h"
#include "Properties/CreateProperties.h"
#include "Properties/Properties.h"

namespace
{
std::unique_ptr<MaterialPropertyLib::Property> createProperty(
    int const geometry_dimension,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    using namespace MaterialPropertyLib;
    // Parsing the property type:
    //! \ogs_file_param{properties__property__type}
    auto const property_type = config.peekConfigParameter<std::string>("type");

    if (property_type == "Constant")
    {
        return createConstant(config);
    }
    if (property_type == "Curve")
    {
        return createCurve(config, curves);
    }
    if (property_type == "Linear")
    {
        return createLinear(config);
    }

    if (property_type == "Exponential")
    {
        return createExponential(config);
    }
    if (property_type == "Function")
    {
        return createFunction(config);
    }

    if (property_type == "Parameter")
    {
        return createParameterProperty(config, parameters);
    }

    if (boost::iequals(property_type, "AverageMolarMass"))
    {
        return createAverageMolarMass(config);
    }

    if (boost::iequals(property_type, "ClausiusClapeyron"))
    {
        return createClausiusClapeyron(config);
    }

    if (boost::iequals(property_type, "CubicLawPermeability"))
    {
        return createCubicLawPermeability(config, parameters);
    }

    if (boost::iequals(property_type, "Dupuit"))
    {
        return createDupuitPermeability(config, parameters);
    }

    if (boost::iequals(property_type,
                       "EffectiveThermalConductivityPorosityMixing"))
    {
        return createEffectiveThermalConductivityPorosityMixing(
            geometry_dimension, config, local_coordinate_system);
    }

    if (boost::iequals(property_type, "IdealGasLaw"))
    {
        return createIdealGasLaw(config);
    }

    if (boost::iequals(property_type, "IdealGasLawBinaryMixture"))
    {
        return createIdealGasLawBinaryMixture(config);
    }

    if (boost::iequals(property_type, "StrainDependentPermeability"))
    {
        return createStrainDependentPermeability(
            geometry_dimension, config, parameters, local_coordinate_system);
    }

    if (boost::iequals(property_type, "GasPressureDependentPermeability"))
    {
        return createGasPressureDependentPermeability(
            geometry_dimension, config, parameters, local_coordinate_system);
    }

    if (boost::iequals(property_type, "EmbeddedFracturePermeability"))
    {
        return createEmbeddedFracturePermeability(geometry_dimension, config,
                                                  parameters);
    }

    if (boost::iequals(property_type,
                       "OrthotropicEmbeddedFracturePermeability"))
    {
        return createOrthotropicEmbeddedFracturePermeability(
            geometry_dimension, config, parameters);
    }

    if (boost::iequals(property_type,
                       "PermeabilityMohrCoulombFailureIndexModel"))
    {
        return createPermeabilityMohrCoulombFailureIndexModel(
            geometry_dimension, config, parameters, local_coordinate_system);
    }

    if (boost::iequals(property_type, "KozenyCarman"))
    {
        return createKozenyCarmanModel(config, parameters);
    }

    if (boost::iequals(property_type, "VermaPruess"))
    {
        return createVermaPruessModel(config, parameters);
    }

    if (boost::iequals(property_type, "PermeabilityOrthotropicPowerLaw"))
    {
        return createPermeabilityOrthotropicPowerLaw(config,
                                                     local_coordinate_system);
    }

    if (boost::iequals(property_type, "PorosityFromMassBalance"))
    {
        return createPorosityFromMassBalance(config, parameters);
    }

    if (boost::iequals(property_type, "TransportPorosityFromMassBalance"))
    {
        return createTransportPorosityFromMassBalance(config, parameters);
    }

    if (boost::iequals(property_type, "SaturationBrooksCorey"))
    {
        return createSaturationBrooksCorey(config);
    }

    if (boost::iequals(property_type, "RelPermBrooksCorey"))
    {
        return createRelPermBrooksCorey(config);
    }

    if (boost::iequals(property_type, "RelPermBrooksCoreyNonwettingPhase"))
    {
        return createRelPermBrooksCoreyNonwettingPhase(config);
    }

    if (boost::iequals(property_type, "SaturationLiakopoulos"))
    {
        return createSaturationLiakopoulos(config);
    }

    if (boost::iequals(property_type, "RelPermLiakopoulos"))
    {
        return createRelPermLiakopoulos(config);
    }

    if (boost::iequals(property_type, "SaturationExponential"))
    {
        return createSaturationExponential(config);
    }
    if (boost::iequals(property_type, "SaturationVanGenuchten"))
    {
        return createSaturationVanGenuchten(config);
    }

    if (boost::iequals(property_type, "CapillaryPressureVanGenuchten"))
    {
        return createCapillaryPressureVanGenuchten(config);
    }

    if (boost::iequals(property_type,
                       "CapillaryPressureRegularizedVanGenuchten"))
    {
        return createCapillaryPressureRegularizedVanGenuchten(config);
    }

    if (boost::iequals(property_type, "RelativePermeabilityVanGenuchten"))
    {
        return createRelPermVanGenuchten(config);
    }

    if (boost::iequals(property_type,
                       "RelativePermeabilityNonWettingPhaseVanGenuchtenMualem"))
    {
        return createRelPermNonWettingPhaseVanGenuchtenMualem(config);
    }

    if (boost::iequals(property_type, "RelativePermeabilityUdell"))
    {
        return createRelPermUdell(config);
    }

    if (boost::iequals(property_type,
                       "RelativePermeabilityUdellNonwettingPhase"))
    {
        return createRelPermUdellNonwettingPhase(config);
    }

    if (boost::iequals(property_type, "SaturationDependentSwelling"))
    {
        return createSaturationDependentSwelling(config,
                                                 local_coordinate_system);
    }

    if (boost::iequals(property_type, "TemperatureDependentFraction"))
    {
        return createTemperatureDependentFraction(config);
    }

    if (boost::iequals(property_type, "SpecificHeatCapacityWithLatentHeat"))
    {
        return createSpecificHeatCapacityWithLatentHeat(config);
    }

    if (boost::iequals(property_type, "BishopsPowerLaw"))
    {
        return createBishopsPowerLaw(config);
    }

    if (boost::iequals(property_type, "BishopsSaturationCutoff"))
    {
        return createBishopsSaturationCutoff(config);
    }

    if (boost::iequals(property_type, "LinearSaturationSwellingStress"))
    {
        return createLinearSaturationSwellingStress(config);
    }

    if (boost::iequals(property_type, "SaturationWeightedThermalConductivity"))
    {
        return createSaturationWeightedThermalConductivity(geometry_dimension,
                                                           config, parameters);
    }

    if (boost::iequals(property_type, "WaterDensityIAPWSIF97Region1"))
    {
        return createWaterDensityIAPWSIF97Region1(config);
    }

    if (boost::iequals(property_type,
                       "WaterSaturationTemperatureIAPWSIF97Region4"))
    {
        return createWaterSaturationTemperatureIAPWSIF97Region4(config);
    }

    if (boost::iequals(property_type, "WaterEnthalpyIAPWSIF97Region1"))
    {
        return createWaterEnthalpyIAPWSIF97Region1(config);
    }

    if (boost::iequals(property_type, "WaterLiquidDensityIAPWSIF97Region4"))
    {
        return createWaterLiquidDensityIAPWSIF97Region4(config);
    }

    if (boost::iequals(property_type, "WaterVapourDensityIAPWSIF97Region4"))
    {
        return createWaterVapourDensityIAPWSIF97Region4(config);
    }

    if (boost::iequals(property_type, "WaterLiquidEnthalpyIAPWSIF97Region4"))
    {
        return createWaterLiquidEnthalpyIAPWSIF97Region4(config);
    }

    if (boost::iequals(property_type, "WaterVapourEnthalpyIAPWSIF97Region4"))
    {
        return createWaterVapourEnthalpyIAPWSIF97Region4(config);
    }

    if (boost::iequals(property_type, "WaterTemperatureIAPWSIF97Region1"))
    {
        return createWaterTemperatureIAPWSIF97Region1(config);
    }

    if (boost::iequals(property_type, "WaterVapourDensity"))
    {
        return createWaterVapourDensity(config);
    }

    if (boost::iequals(property_type, "VapourDiffusionFEBEX"))
    {
        return createVapourDiffusionFEBEX(config);
    }

    if (boost::iequals(property_type, "VapourDiffusionPMQ"))
    {
        return createVapourDiffusionPMQ(config);
    }

    if (boost::iequals(property_type, "LinearWaterVapourLatentHeat"))
    {
        return createLinearWaterVapourLatentHeat(config);
    }

    if (boost::iequals(property_type,
                       "WaterVapourLatentHeatWithCriticalTemperature"))
    {
        return createWaterVapourLatentHeatWithCriticalTemperature(config);
    }

    if (boost::iequals(property_type, "TemperatureDependentDiffusion"))
    {
        return createTemperatureDependentDiffusion(config, parameters);
    }

    if (boost::iequals(property_type, "VolumeFractionAverage"))
    {
        return createVolumeFractionAverage(config);
    }

    if (boost::iequals(property_type, "WaterViscosityIAPWS"))
    {
        return createWaterViscosityIAPWS(config);
    }

    if (boost::iequals(property_type, "LiquidViscosityVogels"))
    {
        return createLiquidViscosityVogels(config);
    }

    if (boost::iequals(property_type, "WaterThermalConductivityIAPWS"))
    {
        return createWaterThermalConductivityIAPWS(config);
    }

    // If none of the above property types are found, OGS throws an error.
    OGS_FATAL("The specified component property type '{:s}' was not recognized",
              property_type);
}
}  // namespace

namespace MaterialPropertyLib
{
std::unique_ptr<PropertyArray> createProperties(
    int const geometry_dimension,
    std::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    if (!config)
    {
        return nullptr;
    }

    //! \ogs_file_param{properties__property}
    auto const& property_configs = config->getConfigSubtreeList("property");
    if (property_configs.empty())
    {
        return nullptr;
    }

    auto properties = std::make_unique<PropertyArray>();

    for (auto property_config : property_configs)
    {
        // Parsing the property name:
        auto const property_name =
            //! \ogs_file_param{properties__property__name}
            property_config.getConfigParameter<std::string>("name");
        // Create a new property based on the configuration subtree:
        auto property =
            createProperty(geometry_dimension, property_config, parameters,
                           local_coordinate_system, curves);

        // Insert the new property at the right position into the components
        // private PropertyArray:
        (*properties)[convertStringToProperty(property_name)] =
            std::move(property);
    }
    return properties;
}

}  // namespace MaterialPropertyLib
