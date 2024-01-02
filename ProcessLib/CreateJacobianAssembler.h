/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <optional>

#include "BaseLib/ConfigTree-fwd.h"

namespace ProcessLib
{
class AbstractJacobianAssembler;

std::unique_ptr<AbstractJacobianAssembler> createJacobianAssembler(
    std::optional<BaseLib::ConfigTree> const& config);
}  // namespace ProcessLib
