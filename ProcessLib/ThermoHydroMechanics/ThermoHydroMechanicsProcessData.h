/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>
#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"
#include "MaterialLib/Fluid/FluidProperty.h"

namespace MeshLib
{
class Element;
}

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
template <typename T>
struct Parameter;

namespace ThermoHydroMechanics
{
template <int DisplacementDim>
struct ThermoHydroMechanicsProcessData
{
    using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

    ThermoHydroMechanicsProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material_,
        std::unique_ptr<MaterialLib::Fluid::FluidProperties>&&
            fluid_properties_,
        Parameter<double> const& intrinsic_permeability_,
        Parameter<double> const& specific_storage_,
        Parameter<double> const& biot_coefficient_,
        Parameter<double> const& porosity_,
        Parameter<double> const& solid_density_,
        Parameter<double> const& solid_linear_thermal_expansion_coefficient_,
        Parameter<double> const&
            fluid_volumetric_thermal_expansion_coefficient_,
        Parameter<double> const& solid_specific_heat_capacity_,
        Parameter<double> const& solid_thermal_conductivity_,
        Parameter<double> const& reference_temperature_,
        Eigen::Matrix<double, DisplacementDim, 1>
            specific_body_force_)
        : material{std::move(material_)},
          fluid_properties{std::move(fluid_properties_)},
          intrinsic_permeability(intrinsic_permeability_),
          specific_storage(specific_storage_),
          biot_coefficient(biot_coefficient_),
          porosity(porosity_),
          solid_density(solid_density_),
          solid_linear_thermal_expansion_coefficient(
              solid_linear_thermal_expansion_coefficient_),
          fluid_volumetric_thermal_expansion_coefficient(
              fluid_volumetric_thermal_expansion_coefficient_),
          solid_specific_heat_capacity(solid_specific_heat_capacity_),
          solid_thermal_conductivity(solid_thermal_conductivity_),
          reference_temperature(reference_temperature_),
          specific_body_force(std::move(specific_body_force_))
    {
    }

    ThermoHydroMechanicsProcessData(ThermoHydroMechanicsProcessData&& other)
        : material{std::move(other.material)},
          fluid_properties{std::move(other.fluid_properties)},
          intrinsic_permeability(other.intrinsic_permeability),
          specific_storage(other.specific_storage),
          biot_coefficient(other.biot_coefficient),
          porosity(other.porosity),
          solid_density(other.solid_density),
          solid_linear_thermal_expansion_coefficient(
              other.solid_linear_thermal_expansion_coefficient),
          fluid_volumetric_thermal_expansion_coefficient(
              other.fluid_volumetric_thermal_expansion_coefficient),
          solid_specific_heat_capacity(other.solid_specific_heat_capacity),
          solid_thermal_conductivity(other.solid_thermal_conductivity),
          reference_temperature(other.reference_temperature),
          specific_body_force(other.specific_body_force),
          dt(other.dt),
          t(other.t)
    {
    }

    //! Copies are forbidden.
    ThermoHydroMechanicsProcessData(ThermoHydroMechanicsProcessData const&) =
        delete;

    //! Assignments are not needed.
    void operator=(ThermoHydroMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(ThermoHydroMechanicsProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    std::unique_ptr<MaterialLib::Fluid::FluidProperties> fluid_properties;
    Parameter<double> const& intrinsic_permeability;
    Parameter<double> const& specific_storage;
    Parameter<double> const& biot_coefficient;
    Parameter<double> const& porosity;
    Parameter<double> const& solid_density;
    Parameter<double> const& solid_linear_thermal_expansion_coefficient;
    Parameter<double> const& fluid_volumetric_thermal_expansion_coefficient;
    Parameter<double> const& solid_specific_heat_capacity;
    Parameter<double> const& solid_thermal_conductivity;
    Parameter<double> const& reference_temperature;
    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt = 0.0;
    double t = 0.0;
    double getFluidDensity(const double p, const double T) const
    {
        ArrayType vars;
        vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
        vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
        return fluid_properties->getValue(
            MaterialLib::Fluid::FluidPropertyType::Density, vars);
    }

    double getFluidViscosity(const double p, const double T) const
    {
        ArrayType vars;
        vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
        vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
        return fluid_properties->getValue(
            MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);
    }
    double getFluidThermalConductivity(const double p, const double T) const
    {
        ArrayType vars;
        vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
        vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
        return fluid_properties->getValue(
            MaterialLib::Fluid::FluidPropertyType::ThermalConductivity, vars);
    }
    double getFluidHeatCapacity(const double p, const double T) const
    {
        ArrayType vars;
        vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
        vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
        return fluid_properties->getValue(
            MaterialLib::Fluid::FluidPropertyType::HeatCapacity, vars);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
