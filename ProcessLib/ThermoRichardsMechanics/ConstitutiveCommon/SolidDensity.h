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
#include "Porosity.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct SolidDensityData
{
    double rho_SR;
    double dry_density_solid;

    static auto reflect()
    {
        return ProcessLib::Reflection::reflectWithName(
            "dry_density_solid", &SolidDensityData::dry_density_solid);
    }
};

template <int DisplacementDim>
struct SolidDensityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              PorosityData const& poro_data,
              TemperatureData<DisplacementDim> const& T_data,
              SolidDensityData& out) const;
};

extern template struct SolidDensityModel<2>;
extern template struct SolidDensityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
