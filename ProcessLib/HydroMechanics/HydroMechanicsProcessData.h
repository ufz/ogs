// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <memory>
#include <utility>
#include <variant>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/VariableType.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Common/HydroMechanics/CouplingScheme.h"
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

namespace HydroMechanics
{
template <int DisplacementDim>
struct HydroMechanicsProcessData
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

    /// If set mass lumping will be applied to the pressure equation.
    bool const mass_lumping = false;

    Common::HydroMechanics::CouplingScheme coupling_scheme;

    /// ID of hydraulic process.
    int const hydraulic_process_id = 0;

    /// ID of the processes that contains mechanical process.
    int const mechanics_related_process_id = 0;

    const bool use_taylor_hood_elements = false;

    MaterialPropertyLib::Variable const phase_variable{};

    MeshLib::PropertyVector<double>* pressure_interpolated = nullptr;
    std::array<MeshLib::PropertyVector<double>*, 3> principal_stress_vector = {
        nullptr, nullptr, nullptr};
    MeshLib::PropertyVector<double>* principal_stress_values = nullptr;

    /// Total permeability as a symmetric tensor of length 4 or 6
    /// with elements in the order k_xx, k_yy, k_zz, k_xy, k_yz, k_xz
    MeshLib::PropertyVector<double>* permeability = nullptr;

    bool isMonolithicSchemeUsed() const
    {
        return std::holds_alternative<Common::HydroMechanics::Monolithic>(
            coupling_scheme);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace HydroMechanics
}  // namespace ProcessLib
