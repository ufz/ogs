/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace BaseLib
{
class ConfigTree;
}
namespace ProcessLib
{
class SecondaryVariableCollection;
}

namespace ProcessLib
{
void createSecondaryVariables(BaseLib::ConfigTree const& config,
                              SecondaryVariableCollection& secondary_variables);

}  // namespace ProcessLib
