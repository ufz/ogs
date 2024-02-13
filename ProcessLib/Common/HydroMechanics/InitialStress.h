/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 13, 2024, 10:13 AM
 */

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
