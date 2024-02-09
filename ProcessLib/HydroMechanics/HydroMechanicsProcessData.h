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
#include <variant>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/VariableType.h"
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
namespace HydroMechanics
{

struct Monolithic
{
};

struct Staggered
{
    /// An optional input to set an algorithmic parameter of the staggered
    /// scheme. In the HydroMechanics process the fixed-stress split has been
    /// implemented as staggered scheme with a stabilization parameter to be
    /// set. For more details see [user guide -
    /// conventions](https://www.opengeosys.org/docs/userguide/basics/conventions/#fixed-stress-split-for-hydro-mechanical-processes).
    double const fixed_stress_stabilization_parameter;

    /// If set, the volume stress rate is fixed over time step, e.g.
    /// \f[\dot{\sigma}_v^{n} = \dot{\sigma}_v^{n-1}, \f]
    /// where \f$n\f$ is the time step index.
    /// Otherwise, the volume stress rate is fixed over coupling iteration, e.g.
    /// \f[\dot{\sigma}_v^{n, k} = \dot{\sigma}_v^{n, k-1}, \f]
    /// where \f$k\f$ is the coupling iteration index, and
    /// \f[ \dot{()}^{n, k} = \left(()^{n, k} - ()^{n-1}\right)/dt. \f]
    bool const fixed_stress_over_time_step;
};

using CouplingScheme = std::variant<Monolithic, Staggered>;

template <int DisplacementDim>
struct HydroMechanicsProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    /// The constitutive relation for the mechanical part.
    std::map<int, std::unique_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    InitialStress const initial_stress;

    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    /// If set mass lumping will be applied to the pressure equation.
    bool const mass_lumping;

    CouplingScheme coupling_scheme;

    /// ID of hydraulic process.
    int const hydraulic_process_id;

    /// ID of the processes that contains mechanical process.
    int const mechanics_related_process_id;

    const bool use_taylor_hood_elements;

    MaterialPropertyLib::Variable phase_pressure;

    MeshLib::PropertyVector<double>* pressure_interpolated = nullptr;
    std::array<MeshLib::PropertyVector<double>*, 3> principal_stress_vector = {
        nullptr, nullptr, nullptr};
    MeshLib::PropertyVector<double>* principal_stress_values = nullptr;

    /// Total permeability as a symmetric tensor of length 4 or 6
    /// with elements in the order k_xx, k_yy, k_zz, k_xy, k_yz, k_xz
    MeshLib::PropertyVector<double>* permeability = nullptr;

    bool isMonolithicSchemeUsed() const
    {
        return std::holds_alternative<Monolithic>(coupling_scheme);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace HydroMechanics
}  // namespace ProcessLib
