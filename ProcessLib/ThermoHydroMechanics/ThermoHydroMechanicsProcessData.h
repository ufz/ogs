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
#include <memory>
#include <utility>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "NumLib/NumericalStability/NumericalStabilization.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Common/HydroMechanics/InitialStress.h"

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
namespace ThermoHydroMechanics
{
template <int DisplacementDim>
struct ThermoHydroMechanicsProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    /// The constitutive relation for the mechanical part.
    std::map<int, std::unique_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    /// The constitutive relation for the mechanical part of ice. When no
    /// freezing is configured, this is a nullptr.
    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        ice_constitutive_relation;

    InitialStress const initial_stress;

    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    NumLib::NumericalStabilization stabilizer;

    MeshLib::PropertyVector<double>* element_fluid_density = nullptr;
    MeshLib::PropertyVector<double>* element_viscosity = nullptr;
    MeshLib::PropertyVector<double>* element_stresses = nullptr;
    MeshLib::PropertyVector<double>* pressure_interpolated = nullptr;
    MeshLib::PropertyVector<double>* temperature_interpolated = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
