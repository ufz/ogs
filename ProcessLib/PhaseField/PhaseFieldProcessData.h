// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <memory>
#include <utility>

#include "MaterialLib/SolidModels/PhaseFieldBase.h"
#include "MeshLib/PropertyVector.h"
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
template <typename T>
struct Parameter;

namespace PhaseField
{
template <int DisplacementDim>
struct PhaseFieldProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    std::map<
        int,
        std::shared_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;
    ParameterLib::Parameter<double> const& residual_stiffness;
    ParameterLib::Parameter<double> const& crack_resistance;
    ParameterLib::Parameter<double> const& crack_length_scale;
    ParameterLib::Parameter<double> const& solid_density;
    /// Optional, initial stress field. A symmetric tensor, short vector
    /// representation of length 4 or 6, ParameterLib::Parameter<double>.
    InitialStress const initial_stress;

    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    bool pressurized_crack = false;
    bool propagating_pressurized_crack = false;
    bool static_pressurized_crack = false;
    double irreversible_threshold;
    MaterialLib::Solids::Phasefield::PhaseFieldModel phasefield_model;
    MaterialLib::Solids::Phasefield::EnergySplitModel energy_split_model;
    MaterialLib::Solids::Phasefield::SofteningCurve softening_curve;
    double characteristic_length;
    std::unique_ptr<MaterialLib::Solids::Phasefield::DegradationDerivative>
        degradation_derivative;

    double const unity_pressure = 1.0;
    double pressure = 0.0;
    double pressure_old = 0.0;
    double pressure_error = 0.0;
    double injected_volume = 0.0;
    double crack_volume = 0.0;
    double elastic_energy = 0.0;
    double surface_energy = 0.0;
    double pressure_work = 0.0;
};

}  // namespace PhaseField
}  // namespace ProcessLib
