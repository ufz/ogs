/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"

namespace ProcessLib::SmallDeformation
{
struct SolidDensityData
{
    double rho_SR;
};

template <int DisplacementDim>
struct SolidDensityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              TemperatureData<DisplacementDim> const& T_data,
              SolidDensityData& out) const;
};

extern template struct SolidDensityModel<2>;
extern template struct SolidDensityModel<3>;
}  // namespace ProcessLib::SmallDeformation
