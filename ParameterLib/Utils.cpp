/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Utils.h"

namespace ParameterLib
{
ParameterBase* findParameterByName(
    std::string const& parameter_name,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    // Find corresponding parameter by name.
    auto const it = std::find_if(
        parameters.cbegin(), parameters.cend(),
        [&parameter_name](std::unique_ptr<ParameterBase> const& p) {
            return p->name == parameter_name;
        });

    if (it == parameters.end())
    {
        return nullptr;
    }

    DBUG("Found parameter `{:s}'.", (*it)->name);
    return it->get();
}
}  // namespace ParameterLib
