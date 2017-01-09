/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_NAMEDFUNCTIONPROVIDER_H
#define NUMLIB_NAMEDFUNCTIONPROVIDER_H

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

#endif // NUMLIB_NAMEDFUNCTIONPROVIDER_H
