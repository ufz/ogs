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
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenVector.h"

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
struct SimplifiedElasticityModel
{
    virtual ~SimplifiedElasticityModel() = default;
    virtual double storageContribution(
        MaterialPropertyLib::Phase const&,
        MaterialPropertyLib::VariableArray const&,
        ParameterLib::SpatialPosition const&, double const, double const) = 0;
    virtual double thermalExpansivityContribution(
        Eigen::Matrix<double, 3, 3> const&, MaterialPropertyLib::Phase const&,
        MaterialPropertyLib::VariableArray const&,
        ParameterLib::SpatialPosition const&, double const, double const) = 0;

    static inline auto bulkCompressibilityFromYoungsModulus(
        MaterialPropertyLib::Phase const& solid_phase,
        MaterialPropertyLib::VariableArray const& variables,
        ParameterLib::SpatialPosition const& x_position, double const t,
        double const dt)
    {
        // assuming: nu[0]=nu(1,2), nu[1]=nu(2,3), nu[2]=nu(1,3)
        if (!solid_phase.hasProperty(
                MaterialPropertyLib::PropertyType::youngs_modulus))
        {
            return 0.0;
        }
        auto const E = MaterialPropertyLib::formEigenVector<3>(
            solid_phase[MaterialPropertyLib::PropertyType::youngs_modulus]
                .value(variables, x_position, t, dt));
        auto const nu = MaterialPropertyLib::formEigenVector<3>(
            solid_phase[MaterialPropertyLib::PropertyType::poissons_ratio]
                .value(variables, x_position, t, dt));
        return (E[0] * E[1] + E[0] * E[2] * (1 - 2 * nu[1]) +
                E[1] * E[2] * (1 - 2 * nu[0] - 2 * nu[2])) /
               (E[0] * E[1] * E[2]);
    }
};

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
