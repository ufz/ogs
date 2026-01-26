// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "BishopsData.h"
#include "CapillaryPressureData.h"
#include "MediaData.h"
#include "PorosityData.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "SaturationData.h"
#include "SolidCompressibilityData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
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
