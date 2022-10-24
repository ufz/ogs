/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct BiotData
{
    double alpha;
};

struct BiotModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              BiotData& out) const;
};
}  // namespace ProcessLib::ThermoRichardsMechanics
