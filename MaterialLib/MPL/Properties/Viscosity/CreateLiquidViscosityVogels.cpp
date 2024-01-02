/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 4, 2021, 4:38 PM
 */

#include "CreateLiquidViscosityVogels.h"

#include "BaseLib/ConfigTree.h"
#include "LiquidViscosityVogels.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createLiquidViscosityVogels(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "LiquidViscosityVogels");
    INFO("Using Vogels model, which gives viscosity in SI unit, Pa s");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    auto const fluid_type =
        //! \ogs_file_param{properties__property__LiquidViscosityVogels__liquid_type}
        config.peekConfigParameter<std::string>("liquid_type");

    if (fluid_type == "Water")
    {
        //! \ogs_file_param{properties__property__LiquidViscosityVogels__liquid_type}
        config.checkConfigParameter("liquid_type", "Water");

        const VogelsViscosityConstantsWater constants;
        return std::make_unique<
            LiquidViscosityVogels<VogelsViscosityConstantsWater>>(
            std::move(property_name), std::move(constants));
    }
    if (fluid_type == "CO2")
    {
        //! \ogs_file_param{properties__property__LiquidViscosityVogels__liquid_type}
        config.checkConfigParameter("liquid_type", "CO2");
        const VogelsViscosityConstantsCO2 constants;
        return std::make_unique<
            LiquidViscosityVogels<VogelsViscosityConstantsCO2>>(
            std::move(property_name), std::move(constants));
    }
    if (fluid_type == "CH4")
    {
        //! \ogs_file_param{properties__property__LiquidViscosityVogels__liquid_type}
        config.checkConfigParameter("liquid_type", "CH4");
        const VogelsViscosityConstantsCH4 constants;
        return std::make_unique<
            LiquidViscosityVogels<VogelsViscosityConstantsCH4>>(
            std::move(property_name), std::move(constants));
    }
    OGS_FATAL(
        "The fluid type {:s} for Vogels model is unavailable.\n"
        "The available fluid types are Water, CO2 and CH4\n",
        fluid_type.data());
}
}  // namespace MaterialPropertyLib
