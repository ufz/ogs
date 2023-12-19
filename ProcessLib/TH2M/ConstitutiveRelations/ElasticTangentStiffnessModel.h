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
#include "ElasticTangentStiffnessData.h"
#include "SolidMechanics.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
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
              TemperatureData const& T_data,
              ElasticTangentStiffnessData<DisplacementDim>& out) const;

private:
    SolidConstitutiveRelation<DisplacementDim> const& solid_material_;
};

extern template struct ElasticTangentStiffnessModel<2>;
extern template struct ElasticTangentStiffnessModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
