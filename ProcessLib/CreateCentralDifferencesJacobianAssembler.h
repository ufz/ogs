// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
class AbstractJacobianAssembler;
}

namespace ProcessLib
{
std::unique_ptr<AbstractJacobianAssembler>
createCentralDifferencesJacobianAssembler(BaseLib::ConfigTree const& config);
}  // namespace ProcessLib
