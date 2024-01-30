/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
    double beta_T_SR;  /// Isotropic solid phase volumetric thermal expansion
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
