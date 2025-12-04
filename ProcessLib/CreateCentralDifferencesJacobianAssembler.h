/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on 2025-12-04 14:04:41
 */
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
