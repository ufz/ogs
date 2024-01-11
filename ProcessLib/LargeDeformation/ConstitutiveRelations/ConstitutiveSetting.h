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

namespace ProcessLib::LargeDeformation
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct ConstitutiveSetting
{
    using GradientVectorType = Eigen::Matrix<
        double,
        DisplacementDim * DisplacementDim + (DisplacementDim == 2 ? 1 : 0), 1>;

    /// Evaluate the constitutive setting.
    void eval(ConstitutiveModels<DisplacementDim>& models, double const t,
              double const dt, ParameterLib::SpatialPosition const& x_position,
              MaterialPropertyLib::Medium const& medium, double const T_ref,
              DeformationGradientData<DisplacementDim> const&
                  deformation_gradient_data,
              GradientVectorType const& deformation_gradient_prev,
              StatefulData<DisplacementDim>& state,
              StatefulDataPrev<DisplacementDim> const& prev_state,
              MaterialStateData<DisplacementDim>& mat_state,
              ConstitutiveTempData<DisplacementDim>& tmp,
              ConstitutiveData<DisplacementDim>& cd) const;
};

extern template struct ConstitutiveSetting<2>;
extern template struct ConstitutiveSetting<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::LargeDeformation
