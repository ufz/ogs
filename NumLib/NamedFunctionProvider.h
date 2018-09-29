/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NamedFunction.h"

namespace NumLib
{

//! Interface used for providing named functions.
class NamedFunctionProvider
{
public:
    virtual std::vector<NamedFunction> getNamedFunctions() const
    {
        return std::vector<NamedFunction>{};
    }

    virtual ~NamedFunctionProvider() = default;
};

} // namespace NumLib
