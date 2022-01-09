/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
        return solid_phase[MaterialPropertyLib::PropertyType::storage_contribution]
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
