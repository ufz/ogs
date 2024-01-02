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

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace LiquidFlow
{
struct LiquidFlowData final
{
    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    /// A vector of the rotation matrices for all elements.
    std::vector<Eigen::MatrixXd> const element_rotation_matrices;

    int const mesh_space_dimension;

    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;

    /// It stores aperture size, which is the thickness of 2D element or the
    /// cross section area of 1D element. For 3D element, the value is set to 1.
    ParameterLib::Parameter<double> const& aperture_size;
};

}  // namespace LiquidFlow
}  // namespace ProcessLib
