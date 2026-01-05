// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
