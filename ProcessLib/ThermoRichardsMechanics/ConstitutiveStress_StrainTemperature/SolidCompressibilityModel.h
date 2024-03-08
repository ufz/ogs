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
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Biot.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidCompressibilityData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim, typename SolidMaterial>
struct SolidCompressibilityModel
{
    explicit SolidCompressibilityModel(SolidMaterial const& solid_material)
        : solid_material_(solid_material)
    {
    }

    void eval(SpaceTimeData const& x_t,
              BiotData const& biot_data,
              ElasticTangentStiffnessData<DisplacementDim> const& C_el_data,
              SolidCompressibilityData& out) const
    {
        out.beta_SR = (1 - biot_data()) / solid_material_.getBulkModulus(
                                              x_t.t, x_t.x, &C_el_data.C_el);
    }

private:
    SolidMaterial const& solid_material_;
};
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
