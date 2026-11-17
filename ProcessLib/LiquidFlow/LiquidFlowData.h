// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/VariableType.h"
#include "NumLib/Fem/ShapeMatrixCache.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace LiquidFlow
{
struct LiquidFlowData final
{
    /// This indicates whether the governing equation is a volume balance or a
    /// mass balance. Its value can be `volume` or `mass`. By default, it is
    /// set to `volume`.
    ///
    /// The choice sets the unit of the boundary and surface fluxes:
    ///
    /// - `volume`: the phase density must be constant. The Neumann boundary
    ///   condition, the surface flux of \c calculatesurfaceflux and the flux
    ///   compared against the threshold of a \c ConstraintDirichlet boundary
    ///   condition are the volumetric (Darcy) flux
    ///   \f$\mathbf{q}\cdot\mathbf{n}\f$ in m/s. The nodal residual output
    ///   is named \c VolumetricFlowRate.
    /// - `mass`: the phase density may depend on pressure.
    ///   The Neumann boundary condition, the surface flux and the
    ///   \c ConstraintDirichlet flux are the mass flux
    ///   \f$\rho\,\mathbf{q}\cdot\mathbf{n}\f$ in kg/m^3*m/s = kg/m^2/s.
    ///   The nodal residual output is named \c MassFlowRate.
    ///
    /// The secondary variable Darcy velocity \c v is the volumetric flux in
    /// m/s for both types.
    bool const is_volume_balance_equation_type;

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

    MaterialPropertyLib::Variable const phase_variable;

    /// If true, the storage (mass) matrix is lumped: column sums for linear
    /// elements, HRZ lumping for higher order elements.
    bool const apply_mass_lumping;
};

}  // namespace LiquidFlow
}  // namespace ProcessLib
