/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CreateCapillaryPressureModel.cpp
 *
 * Created on November 1, 2016, 10:06 AM
 */

#include "CreateCapillaryPressureModel.h"

#include <array>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "CapillaryPressureSaturation.h"
#include "BrookCoreyCapillaryPressureSaturation.h"
#include "vanGenuchtenCapillaryPressureSaturation.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
    \param config ConfigTree object which contains the input data
                  including <type>BrookCorey</type>
                  and it has a tag of <capillary_pressure>
*/
static std::unique_ptr<CapillaryPressureSaturation>
                         createBrookCorey(BaseLib::ConfigTree const& config)
{
    std::array<double, 5> parameters = {
        {//! \ogs_file_param{material_property__porous_medium__porous_medium__capillary_pressure__type__pd}
         config.getConfigParameter<double>("pd"),
         //! \ogs_file_param{material_property__porous_medium__porous_medium__capillary_pressure__type__sr}
         config.getConfigParameter<double>("sr"),
         //! \ogs_file_param{material_property__porous_medium__porous_medium__capillary_pressure__type__smax}
         config.getConfigParameter<double>("smax"),
         //! \ogs_file_param{material_property__porous_medium__porous_medium__capillary_pressure__type__m}
         config.getConfigParameter<double>("m"),
         //! \ogs_file_param{material_property__porous_medium__porous_medium__capillary_pressure__type__pc_max}
         config.getConfigParameter<double>("pc_max")}};

    assert(parameters[3] >= 1.0);  // m >= 1
    return std::unique_ptr<CapillaryPressureSaturation>(
                        new BrookCoreyCapillaryPressureSaturation(parameters));
}

/**
    \param config ConfigTree object which contains the input data
                  including <type>vanGenuchten</type>
                  and it has a tag of <capillary_pressure>
*/
static std::unique_ptr<CapillaryPressureSaturation>
                         createVanGenuchten(BaseLib::ConfigTree const& config)
{
    std::array<double, 5> parameters = {
        {//! \ogs_file_param{material_property__porous_medium__porous_medium__capillary_pressure__type__pd}
         config.getConfigParameter<double>("pd"),
         //! \ogs_file_param{material_property__porous_medium__porous_medium__capillary_pressure__type__sr}
         config.getConfigParameter<double>("sr"),
         //! \ogs_file_param{material_property__porous_medium__porous_medium__capillary_pressure__type__smax}
         config.getConfigParameter<double>("smax"),
         //! \ogs_file_param{material_property__porous_medium__porous_medium__capillary_pressure__type__m}
         config.getConfigParameter<double>("m"),
         //! \ogs_file_param{material_property__porous_medium__porous_medium__capillary_pressure__type__pc_max}
         config.getConfigParameter<double>("pc_max")}};

    assert(parameters[3] <= 1.0);  // m <= 1
    return std::unique_ptr<CapillaryPressureSaturation>(
                      new vanGenuchtenCapillaryPressureSaturation(parameters));
}

std::unique_ptr<CapillaryPressureSaturation> createCapillaryPressureModel(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material_property__porous_medium__porous_medium__capillary_pressure__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "BrookCorey")
    {
        return createBrookCorey(config);
    }
    else if (type == "vanGenuchten")
    {
        return createVanGenuchten(config);
    }
    else
    {
        OGS_FATAL(
            "The capillary pressure mode %s is unavailable.\n"
            "The available types are: \n\tBrookCorey, \n\tvanGenuchten.\n",
            type.data());
    }
}

}  // end namespace
}  // end namespace
