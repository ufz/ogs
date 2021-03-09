/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <boost/optional/optional.hpp>

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"
#include "UserPunch.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::unique_ptr<UserPunch> createUserPunch(
    std::optional<BaseLib::ConfigTree> const& config, MeshLib::Mesh const& mesh)
{
    if (!config)
    {
        return nullptr;
    }

    std::vector<SecondaryVariable> secondary_variables;
    for (auto const& variable_name :
         //! \ogs_file_param{prj__chemical_system__user_punch__headline}
         config->getConfigParameter<std::vector<std::string>>("headline"))
    {
        auto value = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh),
            variable_name,
            MeshLib::MeshItemType::IntegrationPoint,
            1);
        std::fill(std::begin(*value),
                  std::end(*value),
                  std::numeric_limits<double>::quiet_NaN());

        secondary_variables.emplace_back(variable_name, value);
    }

    std::vector<std::string> statements;
    for (auto const& statement :
         //! \ogs_file_param{prj__chemical_system__user_punch__statement}
         config->getConfigParameterList<std::string>("statement"))
    {
        statements.emplace_back(statement);
    }

    return std::make_unique<UserPunch>(std::move(secondary_variables),
                                       std::move(statements));
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
