/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Core>
#include <memory>
#include <utility>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ParameterLib/Parameter.h"
#include "PhaseTransitionModels/PhaseTransitionModel.h"
#include "ProcessLib/Common/HydroMechanics/InitialStress.h"

namespace ProcessLib
{
namespace TH2M
{
template <int DisplacementDim>
struct TH2MProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    /// The constitutive relation for the mechanical part.
    std::map<int, std::unique_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    std::unique_ptr<PhaseTransitionModel> phase_transition_model_ = nullptr;

    ParameterLib::Parameter<double> const& reference_temperature;

    InitialStress const initial_stress;

    /// Specific body forces applied to solid and both fluid phases.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    bool const apply_mass_lumping;

    const bool use_TaylorHood_elements;

    MeshLib::PropertyVector<double>* element_saturation = nullptr;
    MeshLib::PropertyVector<double>* gas_pressure_interpolated = nullptr;
    MeshLib::PropertyVector<double>* capillary_pressure_interpolated = nullptr;
    MeshLib::PropertyVector<double>* temperature_interpolated = nullptr;
    MeshLib::PropertyVector<double>* liquid_pressure_interpolated = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace TH2M
}  // namespace ProcessLib
