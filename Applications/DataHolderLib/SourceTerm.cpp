/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

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
        return ConditionType::NODAL;
    else if (str == "Volume")
        return ConditionType::VOLUME;

    return ConditionType::NONE;
}

std::string SourceTerm::convertTypeToString(ConditionType type)
{
    if (type == ConditionType::NODAL)
        return "Nodal";
    else if (type == ConditionType::VOLUME)
        return "Volume";

    return "";
}

}  // namespace DataHolderLib
