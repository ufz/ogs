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

#include <memory>
#include <utility>

#include <Eigen/Eigen>

#include "ParameterLib/Parameter.h"

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
namespace ThermoMechanics
{
template <int DisplacementDim>
struct ThermoMechanicsProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    /// The constitutive relation for the mechanical part.
    std::map<int, std::unique_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    /// Optional, initial stress field. A symmetric tensor, short vector
    /// representation of length 4 or 6, ParameterLib::Parameter<double>.
    ParameterLib::Parameter<double> const* const initial_stress;

    ParameterLib::Parameter<double> const& reference_solid_density;
    ParameterLib::Parameter<double> const& linear_thermal_expansion_coefficient;
    ParameterLib::Parameter<double> const& specific_heat_capacity;
    ParameterLib::Parameter<double> const&
        thermal_conductivity;  // TODO To be changed as a matrix type variable.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    /// ID of the mechanical process.
    int const mechanics_process_id;

    /// ID of heat conduction process.
    int const heat_conduction_process_id;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ThermoMechanics
}  // namespace ProcessLib
