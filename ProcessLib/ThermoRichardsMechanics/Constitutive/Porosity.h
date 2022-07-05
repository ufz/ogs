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

#include "Bishops.h"
#include "MathLib/KelvinVector.h"
#include "Saturation.h"
#include "SolidCompressibility.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct PorosityData
{
    double phi;
};

template <int DisplacementDim>
struct PorosityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              SolidCompressibilityData const& solid_compressibility_data,
              SaturationData const& S_L_data,
              SaturationData const& S_L_prev_data,
              BishopsData const& bishops_data,
              BishopsData const& bishops_data_prev,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              StrainData<DisplacementDim> const& eps_data,
              StrainData<DisplacementDim> const& eps_prev_data,
              PorosityData const& poro_prev_data, PorosityData& out) const;
};

extern template struct PorosityModel<2>;
extern template struct PorosityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
