// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct SolidThermalExpansionData
{
    KelvinVector<DisplacementDim> solid_linear_thermal_expansivity_vector;
    double beta_T_SR = nan;  /// Solid phase volumetric thermal expansion
                             /// coefficient.
};

template <int DisplacementDim>
struct SolidThermalExpansionModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              SolidThermalExpansionData<DisplacementDim>& out) const;
};

extern template struct SolidThermalExpansionModel<2>;
extern template struct SolidThermalExpansionModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
