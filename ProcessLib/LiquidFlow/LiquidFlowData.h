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
#include "NumLib/Fem/ShapeMatrixCache.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace LiquidFlow
{

/// Governing equation balance type
enum class EquationBalanceType
{
    volume,
    mass
};

struct LiquidFlowData final
{
    /// This indicates whether the governing equation is a volume balance or a
    /// mass balance. Its value can be `volume` or `mass`. If it is set to
    /// `volume`, note that the phase density must be constant, and the unit of
    /// the Neumann boundary condition is m/s. Otherwise, the unit of the
    /// Neumann boundary condition is kg/m³·m/s = kg/m²/s. By default, it is set
    /// to `volume`.
    EquationBalanceType const equation_balance_type;

    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    /// A vector of the rotation matrices for all elements.
    std::vector<Eigen::MatrixXd> const element_rotation_matrices;

    int const mesh_space_dimension;

    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;

    /// It stores aperture size, which is the thickness of 2D element or the
    /// cross section area of 1D element. For 3D element, the value is set to 1.
    ParameterLib::Parameter<double> const& aperture_size;

    /// caches for each mesh element type the shape matrix
    NumLib::ShapeMatrixCache shape_matrix_cache;
};

}  // namespace LiquidFlow
}  // namespace ProcessLib
