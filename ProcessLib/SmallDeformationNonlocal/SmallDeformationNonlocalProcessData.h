/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Core>
#include <memory>
#include <utility>

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}  // namespace MaterialLib
namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
template <int DisplacementDim>
struct SmallDeformationNonlocalProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    std::map<int, std::shared_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;
    /// Solid's density. A scalar quantity, ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& solid_density;
    /// Specific body forces applied to the solid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    double const reference_temperature =
        std::numeric_limits<double>::quiet_NaN();
    double const internal_length_squared =
        std::numeric_limits<double>::quiet_NaN();

    double crack_volume_old = 0.0;
    double crack_volume = 0.0;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
