// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
