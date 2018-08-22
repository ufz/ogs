/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Parameter/Parameter.h"

#include <memory>
#include <utility>

#include <Eigen/Dense>

#include "ProcessLib/RichardsFlow/RichardsFlowMaterialProperties.h"

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
namespace RichardsMechanics
{
template <int DisplacementDim>
struct RichardsMechanicsProcessData
{
    RichardsMechanicsProcessData(
        std::unique_ptr<
            ProcessLib::RichardsFlow::RichardsFlowMaterialProperties>&&
            flow_material_,
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            solid_material_,
        Parameter<double> const& intrinsic_permeability_,
        Parameter<double> const& fluid_bulk_modulus_,
        Parameter<double> const& biot_coefficient_,
        Parameter<double> const& solid_density_,
        Parameter<double> const& solid_bulk_modulus_,
        Parameter<double> const& temperature_,
        Eigen::Matrix<double, DisplacementDim, 1>
            specific_body_force_)
        : flow_material{std::move(flow_material_)},
          solid_material{std::move(solid_material_)},
          intrinsic_permeability(intrinsic_permeability_),
          fluid_bulk_modulus(fluid_bulk_modulus_),
          biot_coefficient(biot_coefficient_),
          solid_density(solid_density_),
          solid_bulk_modulus(solid_bulk_modulus_),
          temperature(temperature_),
          specific_body_force(std::move(specific_body_force_))
    {
    }

    RichardsMechanicsProcessData(RichardsMechanicsProcessData&& other) =
        default;

    //! Copies are forbidden.
    RichardsMechanicsProcessData(RichardsMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(RichardsMechanicsProcessData const&) = delete;
    void operator=(RichardsMechanicsProcessData&&) = delete;

    std::unique_ptr<ProcessLib::RichardsFlow::RichardsFlowMaterialProperties>
        flow_material;

    /// The constitutive relation for the mechanical part.
    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        solid_material;
    /// Permeability of the solid. A scalar quantity, Parameter<double>.
    Parameter<double> const& intrinsic_permeability;
    /// Fluid's bulk modulus. A scalar quantity, Parameter<double>.
    Parameter<double> const& fluid_bulk_modulus;
    /// Biot coefficient. A scalar quantity, Parameter<double>.
    Parameter<double> const& biot_coefficient;
    /// Density of the solid. A scalar quantity, Parameter<double>.
    Parameter<double> const& solid_density;
    /// Solid's bulk modulus. A scalar quantity, Parameter<double>.
    Parameter<double> const& solid_bulk_modulus;
    /// Reference temperature for material properties. A scalar quantity,
    /// Parameter<double>.
    Parameter<double> const& temperature;
    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt = 0.0;
    double t = 0.0;

    MeshLib::PropertyVector<double>* element_saturation = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace RichardsMechanics
}  // namespace ProcessLib
