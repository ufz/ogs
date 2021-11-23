/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ParameterLib/Parameter.h"

#include <memory>
#include <utility>

#include <Eigen/Dense>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}
namespace ProcessLib
{
namespace HydroMechanics
{
template <int DisplacementDim>
struct HydroMechanicsProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map = nullptr;

    /// The constitutive relation for the mechanical part.
    /// \note Linear elasticity is the only supported one in the moment.
    std::map<
        int,
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    /// Optional, initial stress field. A symmetric tensor, short vector
    /// representation of length 4 or 6, ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const* const initial_stress;

    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    /// If set mass lumping will be applied to the pressure equation.
    bool const mass_lumping;

    /// An optional input to set an algorithmic parameter of the staggered
    /// scheme. In the HydroMechanics process the fixed-stress split has been
    /// implemented as staggered scheme with a stabilization parameter to be
    /// set. For more details see [user guide -
    /// conventions](https://www.opengeosys.org/docs/userguide/basics/conventions/#fixed-stress-split-for-hydro-mechanical-processes).
    double const fixed_stress_stabilization_parameter;

    /// ID of hydraulic process.
    int const hydraulic_process_id;

    /// ID of the processes that contains mechanical process.
    int const mechanics_related_process_id;

    const bool use_taylor_hood_elements;

    MeshLib::PropertyVector<double>* pressure_interpolated = nullptr;
    std::array<MeshLib::PropertyVector<double>*, 3> principal_stress_vector = {
        nullptr, nullptr, nullptr};
    MeshLib::PropertyVector<double>* principal_stress_values = nullptr;

    /// Total permeability as a symmetric tensor of length 4 or 6
    /// with elements in the order k_xx, k_yy, k_zz, k_xy, k_yz, k_xz
    MeshLib::PropertyVector<double>* permeability = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace HydroMechanics
}  // namespace ProcessLib
