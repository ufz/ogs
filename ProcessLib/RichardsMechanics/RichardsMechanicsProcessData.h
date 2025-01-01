/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Core>
#include <memory>
#include <utility>

#include "ComputeMicroPorosity.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
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
namespace RichardsMechanics
{
template <int DisplacementDim>
struct RichardsMechanicsProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    /// The constitutive relation for the mechanical part.
    std::map<int, std::shared_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    InitialStress const initial_stress;

    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    std::optional<MicroPorosityParameters> micro_porosity_parameters;

    bool const apply_mass_lumping;

    /// If set, improves convergence of the global Newton method in unsaturated
    /// regime.
    bool const explicit_hm_coupling_in_unsaturated_zone;

    MeshLib::PropertyVector<double>* element_saturation = nullptr;
    MeshLib::PropertyVector<double>* element_porosity = nullptr;
    MeshLib::PropertyVector<double>* element_stresses = nullptr;
    MeshLib::PropertyVector<double>* pressure_interpolated = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace RichardsMechanics
}  // namespace ProcessLib
