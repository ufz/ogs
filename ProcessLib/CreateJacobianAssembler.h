/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/ConfigTree.h"

namespace ProcessLib
{
class AbstractJacobianAssembler;

std::unique_ptr<AbstractJacobianAssembler> createJacobianAssembler(
    boost::optional<BaseLib::ConfigTree> const& config);
}  // ProcessLib
