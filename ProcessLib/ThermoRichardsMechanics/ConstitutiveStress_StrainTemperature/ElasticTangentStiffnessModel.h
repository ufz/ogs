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

#include "ElasticTangentStiffnessData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Base.h"
#include "TraitsBase.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim>
struct ElasticTangentStiffnessModel
{
    explicit ElasticTangentStiffnessModel(
        SolidConstitutiveRelation<DisplacementDim> const& solid_material)
        : solid_material_(solid_material)
    {
    }

    void eval(SpaceTimeData const& x_t,
              TemperatureData<DisplacementDim> const& T_data,
              ElasticTangentStiffnessData<DisplacementDim>& out) const;

private:
    SolidConstitutiveRelation<DisplacementDim> const& solid_material_;
};

extern template struct ElasticTangentStiffnessModel<2>;
extern template struct ElasticTangentStiffnessModel<3>;
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
