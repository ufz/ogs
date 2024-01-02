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

#include <Eigen/Core>
#include <utility>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "NumLib/NumericalStability/NumericalStabilization.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace HT
{
struct HTProcessData final
{
    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;
    bool const has_fluid_thermal_expansion;
    ParameterLib::Parameter<double> const& solid_thermal_expansion;
    ParameterLib::Parameter<double> const& biot_constant;

    bool const has_gravity;
    int const heat_transport_process_id;
    int const hydraulic_process_id;

    NumLib::NumericalStabilization stabilizer;

    /// Projected specific body force vector: R * R^T * b.
    std::vector<Eigen::VectorXd> const projected_specific_body_force_vectors;
    int const mesh_space_dimension;

    /// The aperture size is the thickness of 2D element or the
    /// cross section area of 1D element. For 3D element, the value is set to 1.
    ParameterLib::Parameter<double> const& aperture_size =
        ParameterLib::ConstantParameter<double>("constant_one", 1.0);
};

}  // namespace HT
}  // namespace ProcessLib
