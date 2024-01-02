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

#include "Base.h"
#include "BaseLib/StrongType.h"

namespace ProcessLib::ThermoRichardsMechanics
{
using BiotData = BaseLib::StrongType<double, struct BiotTag>;

struct BiotModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              BiotData& out) const;
};
}  // namespace ProcessLib::ThermoRichardsMechanics
