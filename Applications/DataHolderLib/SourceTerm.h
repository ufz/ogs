/**
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

    class SourceTerm : public DataHolderLib::FemCondition
    {
    public:
        SourceTerm(std::string const process_var, std::string const param_name, ConditionType type)
        : FemCondition(process_var, param_name, type)
        {}

    private:

    };


} // namespace
