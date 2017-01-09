/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace BaseLib { class ConfigTree; }
namespace ProcessLib { class SecondaryVariableCollection; }
namespace NumLib { class NamedFunctionCaller; }

namespace ProcessLib
{
void parseSecondaryVariables(
    BaseLib::ConfigTree const& config,
    SecondaryVariableCollection& secondary_variables,
    NumLib::NamedFunctionCaller& named_function_caller);

} // namespace ProcessLib
