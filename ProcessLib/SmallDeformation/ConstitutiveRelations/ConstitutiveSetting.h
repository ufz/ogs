/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "ConstitutiveData.h"
#include "ConstitutiveModels.h"

namespace ProcessLib::SmallDeformation
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct ConstitutiveSetting
{
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
