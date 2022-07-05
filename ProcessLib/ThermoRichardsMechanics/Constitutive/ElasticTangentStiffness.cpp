/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElasticTangentStiffness.h"

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void ElasticTangentStiffnessModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, TemperatureData<DisplacementDim> const& T_data,
    ElasticTangentStiffnessData<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;

    MPL::VariableArray variable_array;
    MPL::VariableArray variable_array_prev;

    auto const null_state = solid_material_.createMaterialStateVariables();

    using KV = KelvinVector<DisplacementDim>;

    variable_array[static_cast<int>(MPL::Variable::stress)].emplace<KV>(
        KV::Zero());
    variable_array[static_cast<int>(MPL::Variable::mechanical_strain)]
        .emplace<KV>(KV::Zero());
    variable_array[static_cast<int>(MPL::Variable::temperature)]
        .emplace<double>(T_data.T);

    variable_array_prev[static_cast<int>(MPL::Variable::stress)].emplace<KV>(
        KV::Zero());
    variable_array_prev[static_cast<int>(MPL::Variable::mechanical_strain)]
        .emplace<KV>(KV::Zero());
    // Compute previous temperature by linearly following T_dot over the past
    // timestep.
    variable_array_prev[static_cast<int>(MPL::Variable::temperature)]
        .emplace<double>(T_data.T - T_data.T_dot * x_t.dt);

    auto&& solution = solid_material_.integrateStress(
        variable_array_prev, variable_array, x_t.t, x_t.x, x_t.dt, *null_state);

    if (!solution)
    {
        OGS_FATAL("Computation of elastic tangent stiffness failed.");
    }

    out.C_el = std::move(std::get<2>(*solution));
}

template struct ElasticTangentStiffnessModel<2>;
template struct ElasticTangentStiffnessModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
