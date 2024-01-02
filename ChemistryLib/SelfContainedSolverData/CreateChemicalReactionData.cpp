/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateChemicalReactionData.h"

#include <boost/algorithm/string/predicate.hpp>

#include "BaseLib/ConfigTree.h"
#include "ChemicalReaction.h"

namespace ChemistryLib
{
namespace SelfContainedSolverData
{
std::vector<std::unique_ptr<ChemicalReaction>> createChemicalReactionData(
    BaseLib::ConfigTree const& config)
{
    std::vector<std::unique_ptr<ChemicalReaction>> chemical_reactions;

    for (
        auto const& reaction_config :
        //! \ogs_file_param{prj__chemical_system__chemical_reactions__chemical_reaction}
        config.getConfigSubtreeList("chemical_reaction"))
    {
        auto const stoichiometric_vector =
            //! \ogs_file_param{prj__chemical_system__chemical_reactions__chemical_reaction__stoichiometric_coefficients}
            reaction_config.getConfigParameter<std::vector<double>>(
                "stoichiometric_coefficients");

        auto const reaction_type =
            //! \ogs_file_param{prj__chemical_system__chemical_reactions__chemical_reaction__reaction_type}
            reaction_config.getConfigParameter<std::string>("reaction_type");
        if (boost::iequals(reaction_type, "FirstOrderReaction"))
        {
            chemical_reactions.emplace_back(std::make_unique<
                                            FirstOrderReaction>(
                stoichiometric_vector,
                //! \ogs_file_param{prj__chemical_system__chemical_reactions__chemical_reaction__first_order_rate_constant}
                reaction_config.getConfigParameter<double>(
                    "first_order_rate_constant")));
        }
    }

    return chemical_reactions;
}
}  // namespace SelfContainedSolverData
}  // namespace ChemistryLib
