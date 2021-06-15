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

#include <Eigen/Dense>
#include <memory>

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
struct SimplifiedElasticityModel;
}
}

namespace MaterialPropertyLib
{
class MaterialSpatialDistributionMap;
}

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
struct ThermoRichardsFlowProcessData
{
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map = nullptr;

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
