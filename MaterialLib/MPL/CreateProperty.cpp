/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateProperty.h"

#include <string>
#include <vector>
#include "BaseLib/ConfigTree.h"

#include "Properties/CreateProperties.h"

#include "Properties/Properties.h"

#include "Component.h"
#include "Medium.h"
#include "Phase.h"

namespace
{
std::unique_ptr<MaterialPropertyLib::Property> createProperty(
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
        return createCurveProperty(config, curves);
    }
    if (property_type == "Linear")
    {
        return createLinearProperty(config);
    }

    if (property_type == "Exponential")
    {
        return createExponentialProperty(config);
    }

    if (property_type == "Parameter")
    {
        return createParameterProperty(config, parameters);
    }

    if (boost::iequals(property_type, "Dupuit"))
    {
        return createDupuitPermeability(config, parameters);
    }

    if (boost::iequals(property_type, "IdealGasLaw"))
    {
        return createIdealGasLaw(config);
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

    if (boost::iequals(property_type, "SaturationLiakopoulos"))
    {
        return createSaturationLiakopoulos(config);
    }

    if (boost::iequals(property_type, "RelPermLiakopoulos"))
    {
        return createRelPermLiakopoulos(config);
    }

    if (boost::iequals(property_type, "SaturationVanGenuchten"))
    {
        return createSaturationVanGenuchten(config);
    }

    if (boost::iequals(property_type, "CapillaryPressureVanGenuchten"))
    {
        return createCapillaryPressureVanGenuchten(config);
    }

    if (boost::iequals(property_type, "RelativePermeabilityVanGenuchten"))
    {
        return createRelPermVanGenuchten(config);
    }

    if (boost::iequals(property_type, "SaturationDependentSwelling"))
    {
        return createSaturationDependentSwelling(config,
                                                 local_coordinate_system);
    }

    if (boost::iequals(property_type, "BishopsPowerLaw"))
    {
        return createBishopsPowerLaw(config);
    }

    if (boost::iequals(property_type, "BishopsSaturationCutoff"))
    {
        return createBishopsSaturationCutoff(config);
    }

    // If none of the above property types are found, OGS throws an error.
    OGS_FATAL("The specified component property type '%s' was not recognized",
              property_type.c_str());
}
}  // namespace

namespace MaterialPropertyLib
{
std::unique_ptr<PropertyArray> createProperties(
    boost::optional<BaseLib::ConfigTree> const& config,
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
        auto property = createProperty(
            property_config, parameters, local_coordinate_system, curves);

        // Insert the new property at the right position into the components
        // private PropertyArray:
        (*properties)[convertStringToProperty(property_name)] =
            std::move(property);
    }
    return properties;
}

}  // namespace MaterialPropertyLib
