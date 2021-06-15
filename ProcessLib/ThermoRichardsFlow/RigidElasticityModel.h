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
struct RigidElasticityModel : SimplifiedElasticityModel
{
    RigidElasticityModel() { DBUG("using rigid mechanics model"); }

    double storageContribution(MaterialPropertyLib::Phase const&,
                               MaterialPropertyLib::VariableArray const&,
                               ParameterLib::SpatialPosition const&,
                               double const, double const) override
    {
        return 0.0;
    }

    double thermalExpansivityContribution(
        Eigen::Matrix<double, 3, 3> const&, MaterialPropertyLib::Phase const&,
        MaterialPropertyLib::VariableArray const&,
        ParameterLib::SpatialPosition const&, double const,
        double const) override
    {
        return 0.0;
    }
};

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
