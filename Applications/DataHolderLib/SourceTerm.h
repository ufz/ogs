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

#pragma once
#include "FemCondition.h"

namespace DataHolderLib
{
/// Managing data associated with a source term
class SourceTerm final : public FemCondition
{
public:
    enum class ConditionType
    {
        NONE,
        NODAL,
        VOLUME
    };

    SourceTerm(ProcessVariable const& process_var,
               std::string const& param_name, ConditionType type);

    std::string const getConditionClassStr() const override
    {
        return "Source Term";
    }

    /// Returns the type of source term this is
    ConditionType getType() const { return _type; }

    /// Converts the type enum into a string
    static std::string convertTypeToString(ConditionType type);

    /// Converts a string specifying the type into an enum
    static ConditionType convertStringToType(std::string const& str);

private:
    ConditionType _type;
};

}  // namespace DataHolderLib
