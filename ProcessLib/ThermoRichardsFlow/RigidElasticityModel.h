// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "SimplifiedElasticityModel.h"

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
struct RigidElasticityModel : SimplifiedElasticityModel
{
    RigidElasticityModel() { DBUG("using rigid mechanics model"); }

    double storageContribution(
        MaterialPropertyLib::Phase const& /*unused*/,
        MaterialPropertyLib::VariableArray const& /*unused*/,
        ParameterLib::SpatialPosition const& /*unused*/,
        double const /*unused*/, double const /*unused*/) override
    {
        return 0.0;
    }

    double thermalExpansivityContribution(
        Eigen::Matrix<double, 3, 3> const& /*unused*/,
        MaterialPropertyLib::Phase const& /*unused*/,
        MaterialPropertyLib::VariableArray const& /*unused*/,
        ParameterLib::SpatialPosition const& /*unused*/,
        double const /*unused*/, double const /*unused*/) override
    {
        return 0.0;
    }
};

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
