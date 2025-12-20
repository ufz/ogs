// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Water.h"

#include "MaterialLib/MPL/Properties/Properties.h"

namespace MaterialPropertyLib
{
Water::Water(std::unique_ptr<PropertyArray>&& properties)
    : Component{"Water", std::move(properties)}
{
}
}  // namespace MaterialPropertyLib
