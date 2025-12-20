// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
