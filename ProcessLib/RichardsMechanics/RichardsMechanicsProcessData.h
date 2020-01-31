/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    std::unique_ptr<ProcessLib::RichardsFlow::RichardsFlowMaterialProperties>
        flow_material;

    /// The constitutive relation for the mechanical part.
    std::map<
        int,
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    /// Optional, initial stress field. A symmetric tensor, short vector
    /// representation of length 4 or 6, ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const* const initial_stress;

    /// Fluid's bulk modulus. A scalar quantity,
    /// ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& fluid_bulk_modulus;
    /// Biot coefficient. A scalar quantity, ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& biot_coefficient;
    /// Density of the solid. A scalar quantity,
    /// ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& solid_density;
    /// Solid's bulk modulus. A scalar quantity,
    /// ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& solid_bulk_modulus;
    /// Reference temperature for material properties. A scalar quantity,
    /// ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const& temperature;
    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    bool const apply_mass_lumping;

    MeshLib::PropertyVector<double>* element_saturation = nullptr;
    MeshLib::PropertyVector<double>* element_stresses = nullptr;
    MeshLib::PropertyVector<double>* pressure_interpolated = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace RichardsMechanics
}  // namespace ProcessLib
