/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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

namespace
{
ChemistryLib::PhreeqcIOData::ChargeBalance parseChargeBalance(
    BaseLib::ConfigTree const& config)
{
    auto const charge_balance_in_str =
        //! \ogs_file_param{prj__chemical_system__solution__charge_balance}
        config.getConfigParameter<std::string>("charge_balance", "");

    if (charge_balance_in_str.empty())
    {
        return ChemistryLib::PhreeqcIOData::ChargeBalance::Unspecified;
    }
    if (charge_balance_in_str == "pH")
    {
        return ChemistryLib::PhreeqcIOData::ChargeBalance::pH;
    }
    if (charge_balance_in_str == "pe")
    {
        return ChemistryLib::PhreeqcIOData::ChargeBalance::pe;
    }

    OGS_FATAL(
        "Error in specifying means of adjusting charge. Achieving charge "
        "balance is currently supported with the way of adjusting pH value or "
        "pe value.");
}
}  // namespace

namespace ChemistryLib
{
namespace PhreeqcIOData
{
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

    auto charge_balance = parseChargeBalance(config);

    return {temperature, pressure, pe, std::move(components), charge_balance};
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
