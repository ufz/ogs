/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateAqueousSolution.h"
#include "AqueousSolution.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "CreateSolutionComponent.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
namespace
{
MeansOfAdjustingCharge parseMeansOfAdjustingCharge(
    BaseLib::ConfigTree const& config)
{
    auto const means_of_adjusting_charge_in_str =
        //! \ogs_file_param{prj__chemical_system__solution__means_of_adjusting_charge}
        config.getConfigParameter<std::string>("means_of_adjusting_charge", "");

    if (means_of_adjusting_charge_in_str.empty())
    {
        return MeansOfAdjustingCharge::Unspecified;
    }
    if (means_of_adjusting_charge_in_str == "pH")
    {
        return MeansOfAdjustingCharge::pH;
    }
    if (means_of_adjusting_charge_in_str == "pe")
    {
        return MeansOfAdjustingCharge::pe;
    }

    OGS_FATAL(
        "Error in specifying means of adjusting charge. Achieving charge "
        "balance is currently supported with the way of adjusting pH value or "
        "pe value.");
}
}  // namespace

AqueousSolution createAqueousSolution(
    BaseLib::ConfigTree const& config,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map)
{
    //! \ogs_file_param{prj__chemical_system__solution__temperature}
    auto const temperature = config.getConfigParameter<double>("temperature");

    //! \ogs_file_param{prj__chemical_system__solution__pressure}
    auto const pressure = config.getConfigParameter<double>("pressure");

    //! \ogs_file_param{prj__chemical_system__solution__pe}
    auto const pe = config.getConfigParameter<double>("pe");

    auto components =
        createSolutionComponents(config, process_id_to_component_name_map);

    auto means_of_adjusting_charge = parseMeansOfAdjustingCharge(config);

    return {temperature, pressure, pe, std::move(components),
            means_of_adjusting_charge};
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
