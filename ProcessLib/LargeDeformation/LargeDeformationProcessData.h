/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Core>
#include <memory>
#include <utility>

#include "ConstitutiveRelations/SolidMechanics.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "ParameterLib/Parameter.h"

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
namespace LargeDeformation
{
template <int DisplacementDim>
struct LargeDeformationProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    std::map<int,
             std::unique_ptr<ConstitutiveRelations::SolidConstitutiveRelation<
                 DisplacementDim>>>
        solid_materials;

    /// Optional, initial stress field. A symmetric tensor, short vector
    /// representation of length 4 or 6, ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const* const initial_stress;

    /// Specific body forces applied to the solid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    ParameterLib::Parameter<double> const* const reference_temperature;

    std::array<MeshLib::PropertyVector<double>*, 3> principal_stress_vector = {
        nullptr, nullptr, nullptr};
    MeshLib::PropertyVector<double>* principal_stress_values = nullptr;
};

}  // namespace LargeDeformation
}  // namespace ProcessLib
