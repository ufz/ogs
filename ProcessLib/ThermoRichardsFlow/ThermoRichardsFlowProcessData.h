// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <memory>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "SimplifiedElasticityModel.h"

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
struct ThermoRichardsFlowProcessData
{
    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of global mesh dimension's length.
    Eigen::VectorXd const specific_body_force;

    bool const apply_mass_lumping;
    std::unique_ptr<SimplifiedElasticityModel> simplified_elasticity = nullptr;

    MeshLib::PropertyVector<double>* element_saturation = nullptr;
    MeshLib::PropertyVector<double>* element_porosity = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
