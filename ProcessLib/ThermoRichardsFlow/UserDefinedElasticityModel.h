// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "SimplifiedElasticityModel.h"

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
struct UserDefinedElasticityModel : SimplifiedElasticityModel
{
    UserDefinedElasticityModel()
    {
        DBUG("using user defined simplified elasticity model");
    }

    double storageContribution(
        MaterialPropertyLib::Phase const& solid_phase,
        MaterialPropertyLib::VariableArray const& variables,
        ParameterLib::SpatialPosition const& x_position, double const t,
        double const dt) override
    {
        return solid_phase
            [MaterialPropertyLib::PropertyType::storage_contribution]
                .template value<double>(variables, x_position, t, dt);
    }
    double thermalExpansivityContribution(
        Eigen::Matrix<double, 3,
                      3> const& /*solid_linear_thermal_expansion_coefficient*/,
        MaterialPropertyLib::Phase const& solid_phase,
        MaterialPropertyLib::VariableArray const& variables,
        ParameterLib::SpatialPosition const& x_position, double const t,
        double const dt) override
    {
        return solid_phase[MaterialPropertyLib::PropertyType::
                               thermal_expansivity_contribution]
            .template value<double>(variables, x_position, t, dt);
    }
};

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
