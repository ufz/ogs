/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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
    MeshLib::PropertyVector<int> const* const material_ids;

    /// The constitutive relation for the mechanical part.
    /// \note Linear elasticity is the only supported one in the moment.
    std::map<
        int,
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;
    /// Permeability of the solid. A scalar quantity,
    /// ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& intrinsic_permeability;
    /// Volumetric average specific storage of the solid and fluid phases.
    /// A scalar quantity, ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& specific_storage;
    /// Fluid's viscosity. A scalar quantity, ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& fluid_viscosity;
    /// Fluid's density. A scalar quantity, ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& fluid_density;
    /// Biot coefficient. A scalar quantity, ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& biot_coefficient;
    /// Porosity of the solid. A scalar quantity,
    /// ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& porosity;
    /// Solid's density. A scalar quantity, ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& solid_density;
    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double const reference_temperature;

    MeshLib::PropertyVector<double>* pressure_interpolated = nullptr;
    double dt = 0.0;
    double t = 0.0;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace HydroMechanics
}  // namespace ProcessLib
