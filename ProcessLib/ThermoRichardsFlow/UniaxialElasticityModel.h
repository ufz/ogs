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
struct UniaxialElasticityModel : SimplifiedElasticityModel
{
    UniaxialElasticityModel()
    {
        DBUG("using uniaxial simplified mechanics model");
    }

    double storageContribution(
        MaterialPropertyLib::Phase const& solid_phase,
        MaterialPropertyLib::VariableArray const& variables,
        ParameterLib::SpatialPosition const& x_position, double const t,
        double const dt) override
    {
        auto const E = MaterialPropertyLib::formEigenVector<3>(
            solid_phase[MaterialPropertyLib::PropertyType::youngs_modulus]
                .value(variables, x_position, t, dt));
        auto const nu = MaterialPropertyLib::formEigenVector<3>(
            solid_phase[MaterialPropertyLib::PropertyType::poissons_ratio]
                .value(variables, x_position, t, dt));
        auto const nu12 = nu[0];
        auto const nu23 = nu[1];
        auto const nu13 = nu[2];
        auto const nu21 = nu12 * E[1] / E[0];
        auto const nu32 = nu23 * E[2] / E[1];
        auto const nu31 = nu13 * E[2] / E[0];
        auto const D = 1 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 -
                       2 * nu12 * nu23 * nu31;
        return D / (E[2] * (1 - nu12 * nu21));
    }

    double thermalExpansivityContribution(
        Eigen::Matrix<double, 3, 3> const& solid_linear_thermal_expansion_coefficient,
        MaterialPropertyLib::Phase const& solid_phase,
        MaterialPropertyLib::VariableArray const& variables,
        ParameterLib::SpatialPosition const& x_position, double const t,
        double const dt) override
    {
        auto const E = MaterialPropertyLib::formEigenVector<3>(
            solid_phase[MaterialPropertyLib::PropertyType::youngs_modulus]
                .value(variables, x_position, t, dt));
        auto const nu = MaterialPropertyLib::formEigenVector<3>(
            solid_phase[MaterialPropertyLib::PropertyType::poissons_ratio]
                .value(variables, x_position, t, dt));
        auto const nu12 = nu[0];
        auto const nu23 = nu[1];
        auto const nu13 = nu[2];
        auto const nu21 = nu12 * E[1] / E[0];
        auto const D = (1 - nu12 * nu21);
        return -(solid_linear_thermal_expansion_coefficient(2, 2) +
                 solid_linear_thermal_expansion_coefficient(0, 0) *
                     (nu13 + nu12 * nu23) / D +
                 solid_linear_thermal_expansion_coefficient(1, 1) *
                     (nu23 + nu13 * nu21) / D);
    }
};

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
