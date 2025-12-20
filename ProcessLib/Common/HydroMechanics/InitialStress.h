// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
struct InitialStress
{
    enum class Type
    {
        Effective,
        Total
    };

    bool isTotalStress() const { return value && (type == Type::Total); }

    ParameterLib::Parameter<double> const* const value = nullptr;
    Type const type = Type::Effective;
};
}  // namespace ProcessLib
