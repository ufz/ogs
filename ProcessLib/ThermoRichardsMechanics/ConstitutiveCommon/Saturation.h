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

namespace ProcessLib::ThermoRichardsMechanics
{
struct SaturationData
{
    double S_L;

    static auto reflect()
    {
        return ProcessLib::Reflection::reflectWithName("saturation",
                                                       &SaturationData::S_L);
    }
};

struct SaturationDataDeriv
{
    double dS_L_dp_cap;
};

template <int DisplacementDim>
struct SaturationModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              SaturationData& S_L_data, SaturationDataDeriv& dS_L_data) const;
};

extern template struct SaturationModel<2>;
extern template struct SaturationModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
