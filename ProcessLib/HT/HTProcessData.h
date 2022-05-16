/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>
#include <memory>
#include <utility>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "NumLib/NumericalStability/NumericalStabilization.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace HT
{
struct HTProcessData final
{
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
    bool const has_fluid_thermal_expansion;
    ParameterLib::Parameter<double> const& solid_thermal_expansion;
    ParameterLib::Parameter<double> const& biot_constant;
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
    int const heat_transport_process_id;
    int const hydraulic_process_id;

    std::unique_ptr<NumLib::NumericalStabilization> stabilizer;
};

}  // namespace HT
}  // namespace ProcessLib
