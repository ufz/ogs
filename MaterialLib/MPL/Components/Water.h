// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/MPL/Component.h"

namespace MaterialPropertyLib
{
/// A class for Water derived from Component.
///
/// This class can holds material constants and default properties of ordinary
/// water.
struct Water final : public Component
{
    explicit Water(std::unique_ptr<PropertyArray>&& properties);
};
}  // namespace MaterialPropertyLib
