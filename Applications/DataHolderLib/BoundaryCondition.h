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
#include "Applications/DataHolderLib/FemCondition.h"

namespace DataHolderLib
{
class BoundaryCondition final : public DataHolderLib::FemCondition
{
public:
    enum ConditionType
    {
        NONE = 0,
        DIRICHLET,
        NONUNIFORMDIRICHLET,
        NEUMANN,
        NONUNIFORMNEUMANN,
        ROBIN
    };

    BoundaryCondition(ProcessVariable const& process_var,
                      std::string const& param_name, ConditionType type);

    ~BoundaryCondition() = default;

    std::string const getConditionClassStr() const
    {
        return "Boundary Condition";
    }

    /// Returns the type of boundary condition this is
    ConditionType getType() const { return _type; }

    /// Converts the type enum into a string
    static ConditionType convertStringToType(std::string const& str);

    /// Converts a string specifying the type into an enum
    static std::string convertTypeToString(ConditionType type);

private:
    ConditionType _type;
};

}  // namespace
