/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
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
createForwardDifferencesJacobianAssembler(BaseLib::ConfigTree const& config);

}  // namespace ProcessLib
