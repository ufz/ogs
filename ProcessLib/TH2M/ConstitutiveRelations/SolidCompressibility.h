/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Biot.h"
#include "ElasticTangentStiffnessData.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
using SolidCompressibilityData =
    BaseLib::StrongType<double, struct SolidCompressibilityDataTag>;

template <int DisplacementDim, typename SolidMaterial>
struct SolidCompressibilityModel
{
    explicit SolidCompressibilityModel(SolidMaterial const& solid_material)
        : solid_material_(solid_material)
    {
    }

    void eval(SpaceTimeData const& x_t,
              BiotData const& biot,
              ElasticTangentStiffnessData<DisplacementDim> const& C_el_data,
              SolidCompressibilityData& out) const
    {
        *out = (1 - biot()) / solid_material_.getBulkModulus(
                                  x_t.t, x_t.x, &C_el_data.stiffness_tensor);
    }

private:
    SolidMaterial const& solid_material_;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
