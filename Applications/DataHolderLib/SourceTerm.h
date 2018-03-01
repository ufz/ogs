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
#include "Applications/DataHolderLib/Condition.h"

namespace DataHolderLib
{

    class SourceTerm : public DataHolderLib::Condition
    {
    public:
        SourceTerm(std::string const process_var, std::string const param_name, ConditionType type)
        : Condition(process_var, param_name, type)
        {}

    private:

    };


} // namespace
