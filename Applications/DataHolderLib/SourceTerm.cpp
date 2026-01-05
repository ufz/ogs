// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SourceTerm.h"

namespace DataHolderLib
{
SourceTerm::SourceTerm(ProcessVariable const& process_var,
                       std::string const& param_name, ConditionType type)
    : FemCondition(process_var, param_name), _type(type)
{
}

SourceTerm::ConditionType SourceTerm::convertStringToType(
    std::string const& str)
{
    if (str == "Nodal")
    {
        return ConditionType::NODAL;
    }
    if (str == "Volume")
    {
        return ConditionType::VOLUME;
    }

    return ConditionType::NONE;
}

std::string SourceTerm::convertTypeToString(ConditionType type)
{
    if (type == ConditionType::NODAL)
    {
        return "Nodal";
    }
    if (type == ConditionType::VOLUME)
    {
        return "Volume";
    }

    return "";
}

}  // namespace DataHolderLib
