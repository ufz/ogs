// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Biot.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidCompressibilityData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidMechanicsDataStateless.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
// TODO remove typename SolidMaterial
template <int DisplacementDim, typename SolidMaterial>
struct SolidCompressibilityModel
{
    explicit SolidCompressibilityModel(SolidMaterial const& solid_material)
        : solid_material_(solid_material)
    {
    }

    void eval(SpaceTimeData const& x_t,
              BiotData const& biot_data,
              SolidMechanicsDataStateless<DisplacementDim> const&
                  s_mech_data_stateless,
              SolidCompressibilityData& out) const
    {
        out.beta_SR =
            (1 - biot_data()) /
            solid_material_.getBulkModulus(
                x_t.t, x_t.x, &s_mech_data_stateless.stiffness_tensor);
    }

    static SolidCompressibilityModel create(SolidMaterial const& solid_material)
    {
        return SolidCompressibilityModel{solid_material};
    }

private:
    SolidMaterial const& solid_material_;
};
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
