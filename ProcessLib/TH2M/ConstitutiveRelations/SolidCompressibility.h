// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
