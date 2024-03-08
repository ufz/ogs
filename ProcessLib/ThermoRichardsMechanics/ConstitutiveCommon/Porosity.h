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

#include "Bishops.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "Saturation.h"
#include "SolidCompressibilityData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct PorosityData
{
    double phi;

    static auto reflect()
    {
        return ProcessLib::Reflection::reflectWithName("porosity",
                                                       &PorosityData::phi);
    }
};

template <int DisplacementDim>
struct PorosityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              SolidCompressibilityData const& solid_compressibility_data,
              SaturationData const& S_L_data,
              PrevState<SaturationData> const& S_L_prev_data,
              BishopsData const& bishops_data,
              PrevState<BishopsData> const& bishops_data_prev,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              StrainData<DisplacementDim> const& eps_data,
              PrevState<StrainData<DisplacementDim>> const& eps_prev_data,
              PrevState<PorosityData> const& poro_prev_data,
              PorosityData& out) const;
};

extern template struct PorosityModel<2>;
extern template struct PorosityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
