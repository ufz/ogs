/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include "BoundaryCondition.h"

namespace DataHolderLib
{
/// Managing data associated with a boundary condition
BoundaryCondition::BoundaryCondition(ProcessVariable const& process_var,
                                     std::string const& param_name,
                                     ConditionType type)
    : FemCondition(process_var, param_name), _type(type)
{
}

BoundaryCondition::ConditionType BoundaryCondition::convertStringToType(
    std::string const& str)
{
    if (str == "Dirichlet")
    {
        return ConditionType::DIRICHLET;
    }
    if (str == "Neumann")
    {
        return ConditionType::NEUMANN;
    }
    if (str == "Robin")
    {
        return ConditionType::ROBIN;
    }

    return ConditionType::NONE;
}

std::string BoundaryCondition::convertTypeToString(ConditionType type)
{
    if (type == ConditionType::DIRICHLET)
    {
        return "Dirichlet";
    }
    if (type == ConditionType::NEUMANN)
    {
        return "Neumann";
    }
    if (type == ConditionType::ROBIN)
    {
        return "Robin";
    }

    return "";
}

}  // namespace DataHolderLib
