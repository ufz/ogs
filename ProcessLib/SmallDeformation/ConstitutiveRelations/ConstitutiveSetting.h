// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ConstitutiveData.h"
#include "ConstitutiveModels.h"
#include "MaterialState.h"

namespace ProcessLib::SmallDeformation
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct ConstitutiveSetting
{
    /// Initialize the constitutive setting and check evaluation order.
    void init();

    /// Evaluate the constitutive setting.
    void eval(ConstitutiveModels<DisplacementDim>& models, double const t,
              double const dt, ParameterLib::SpatialPosition const& x_position,
              MaterialPropertyLib::Medium const& medium, double const T_ref,
              KelvinVector<DisplacementDim> const& eps,
              KelvinVector<DisplacementDim> const& eps_prev,
              StatefulData<DisplacementDim>& state,
              StatefulDataPrev<DisplacementDim> const& prev_state,
              MaterialStateData<DisplacementDim>& mat_state,
              ConstitutiveTempData<DisplacementDim>& tmp,
              OutputData<DisplacementDim>& out,
              ConstitutiveData<DisplacementDim>& cd) const;
};

extern template struct ConstitutiveSetting<2>;
extern template struct ConstitutiveSetting<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::SmallDeformation
