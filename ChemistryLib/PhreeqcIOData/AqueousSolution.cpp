/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AqueousSolution.h"
#include "BaseLib/ConfigTreeUtil.h"

namespace ChemistryLib
{
AqueousSolution createAqueousSolution(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__chemical_system__solution__temperature}
    auto const temperature = config.getConfigParameter<double>("temperature");

    //! \ogs_file_param{prj__chemical_system__solution__pressure}
    auto const pressure = config.getConfigParameter<double>("pressure");

    //! \ogs_file_param{prj__chemical_system__solution__pe}
    auto const pe = config.getConfigParameter<double>("pe");

    //! \ogs_file_param{prj__chemical_system__solution__components}
    auto comp_config = config.getConfigSubtree("components");

    std::vector<Component> components;
    for (
        auto const& component_name :
        //! \ogs_file_param{prj__chemical_system__solution__components__component}
        comp_config.getConfigParameterList<std::string>("component"))
    {
        components.emplace_back(component_name);
    }

    // conversion the variable 'means_of_adjusting_charge' from std::string to
    // enumerate class.
    auto const means_of_adjusting_charge_in_str =
        //! \ogs_file_param{prj__chemical_system__solution__means_of_adjusting_charge}
        config.getConfigParameterOptional<std::string>(
            "means_of_adjusting_charge");

    MeansOfAdjustingCharge means_of_adjusting_charge;
    if (means_of_adjusting_charge_in_str)
    {
        if (*means_of_adjusting_charge_in_str == "pH")
            means_of_adjusting_charge = MeansOfAdjustingCharge::pH;
        else if (*means_of_adjusting_charge_in_str == "pe")
            means_of_adjusting_charge = MeansOfAdjustingCharge::pe;
        else
            OGS_FATAL(
                "Error in specifying means of adjusting charge. Achieving "
                "charge balance is currently supported with the way of "
                "adjusting pH value or pe value.");
    }
    else
    {
        means_of_adjusting_charge = MeansOfAdjustingCharge::Unspecified;
    }

    AqueousSolution aqueous_solution(temperature,
                                     pressure,
                                     pe,
                                     std::move(components),
                                     means_of_adjusting_charge);

    return aqueous_solution;
}
}  // namespace ChemistryLib
