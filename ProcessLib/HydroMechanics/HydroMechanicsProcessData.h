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
#include "MaterialLib/Fluid/FluidType/FluidType.h"

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
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    /// The constitutive relation for the mechanical part.
    /// \note Linear elasticity is the only supported one in the moment.
    std::map<
        int,
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    /// Optional, initial stress field. A symmetric tensor, short vector
    /// representation of length 4 or 6, ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const* const initial_stress;

    /// Permeability of the solid. A scalar quantity,
    /// ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& intrinsic_permeability;
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

    /// Fluid's compressibility. A scalar quantity.
    /// Only used for compressible_fluid fluid_type
    double const fluid_compressibility =
        std::numeric_limits<double>::quiet_NaN();

    /// Reference Temperature. A scalar quantity.
    /// Only used for ideal_gas fluid_type
    double const reference_temperature =
        std::numeric_limits<double>::quiet_NaN();

    /// Specific gas constant. A scalar quantity.
    /// Only used for ideal_gas fluid_type
    double const specific_gas_constant =
        std::numeric_limits<double>::quiet_NaN();

    /// Fluid type. Enumerator with possible values:
    /// incompressible_fluid, compressible_fluid, ideal_gas
    FluidType::Fluid_Type const fluid_type;

    /// will be removed after linking with MPL
    double getFluidDensity(
        double const& t, ParameterLib::SpatialPosition const& x_position,
        double const& p_fr)
    {
        if (fluid_type == FluidType::Fluid_Type::INCOMPRESSIBLE_FLUID ||
            fluid_type == FluidType::Fluid_Type::COMPRESSIBLE_FLUID)
        {
            return fluid_density(t, x_position)[0];
        }
        if (fluid_type == FluidType::Fluid_Type::IDEAL_GAS)
        {
            return p_fr / (specific_gas_constant * reference_temperature);
        }
        OGS_FATAL("unknown fluid type %d", static_cast<int> (fluid_type));
    }

    /// will be removed after linking with MPL
    double getFluidCompressibility(double const& p_fr)
    {
        if (fluid_type == FluidType::Fluid_Type::INCOMPRESSIBLE_FLUID)
        {
            return 0.0;
        }
        if (fluid_type == FluidType::Fluid_Type::COMPRESSIBLE_FLUID)
        {
            return fluid_compressibility;
        }
        if (fluid_type == FluidType::Fluid_Type::IDEAL_GAS)
        {
            return 1.0 / p_fr;
        }
        OGS_FATAL("unknown fluid type %d", static_cast<int> (fluid_type));
    }

    MeshLib::PropertyVector<double>* pressure_interpolated = nullptr;
    std::array<MeshLib::PropertyVector<double>*, 3> principal_stress_vector = {
        nullptr, nullptr, nullptr};
    MeshLib::PropertyVector<double>* principal_stress_values = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace HydroMechanics
}  // namespace ProcessLib
