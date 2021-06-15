/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on August 14, 2020, 10:56 AM
 */

#pragma once

#include "SimplifiedElasticityModel.h"

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
struct HydrostaticElasticityModel : SimplifiedElasticityModel
{
    HydrostaticElasticityModel()
    {
        DBUG("using hydrostatic simplified mechanics model");
    }

    double storageContribution(
        MaterialPropertyLib::Phase const& solid_phase,
        MaterialPropertyLib::VariableArray const& variable_array,
        ParameterLib::SpatialPosition const& pos, double const t,
        double const dt) override
    {
        return bulkCompressibilityFromYoungsModulus(
            solid_phase, variable_array, pos, t, dt);
    }

    double thermalExpansivityContribution(
        Eigen::Matrix<double, 3, 3> const& solid_linear_thermal_expansion_coefficient,
        MaterialPropertyLib::Phase const& /*solid_phase*/,
        MaterialPropertyLib::VariableArray const& /*variable_array*/,
        ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
        double const /*dt*/) override
    {
        return -solid_linear_thermal_expansion_coefficient.trace();
    }
};

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
