/**
 * \file
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
}  // namespace MaterialLib
namespace ProcessLib
{
namespace ThermoHydroMechanics
{
template <int DisplacementDim>
struct ThermoHydroMechanicsProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    /// The constitutive relation for the mechanical part.
    /// \note Linear elasticity is the only supported one in the moment.
    std::map<
        int,
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;
    /// Permeability of the solid. A scalar quantity,
    ///	ParameterLib::Parameter<double>.
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
    ParameterLib::Parameter<double> const&
        solid_linear_thermal_expansion_coefficient;
    ParameterLib::Parameter<double> const&
        fluid_volumetric_thermal_expansion_coefficient;
    ParameterLib::Parameter<double> const& solid_specific_heat_capacity;
    ParameterLib::Parameter<double> const& fluid_specific_heat_capacity;
    ParameterLib::Parameter<double> const& solid_thermal_conductivity;
    ParameterLib::Parameter<double> const& fluid_thermal_conductivity;
    ParameterLib::Parameter<double> const& reference_temperature;
    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    double dt = std::numeric_limits<double>::quiet_NaN();
    double t = std::numeric_limits<double>::quiet_NaN();

    MeshLib::PropertyVector<double>* pressure_interpolated = nullptr;
    MeshLib::PropertyVector<double>* temperature_interpolated = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
