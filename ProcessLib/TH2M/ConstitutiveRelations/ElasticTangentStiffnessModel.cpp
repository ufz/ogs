/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ElasticTangentStiffnessModel.h"

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void ElasticTangentStiffnessModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, TemperatureData const& T_data,
    ElasticTangentStiffnessData<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;

    MPL::VariableArray variable_array;
    MPL::VariableArray variable_array_prev;

    auto const null_state = solid_material_.createMaterialStateVariables();
    solid_material_.initializeInternalStateVariables(x_t.t, x_t.x, *null_state);

    using KV = KelvinVector<DisplacementDim>;

    variable_array.stress.emplace<KV>(KV::Zero());
    variable_array.mechanical_strain.emplace<KV>(KV::Zero());
    variable_array.temperature = T_data.T;

    variable_array_prev.stress.emplace<KV>(KV::Zero());
    variable_array_prev.mechanical_strain.emplace<KV>(KV::Zero());
    variable_array_prev.temperature = T_data.T_prev;

    auto&& solution = solid_material_.integrateStress(
        variable_array_prev, variable_array, x_t.t, x_t.x, x_t.dt, *null_state);

    if (!solution)
    {
        OGS_FATAL("Computation of elastic tangent stiffness failed.");
    }

    out.stiffness_tensor = std::move(std::get<2>(*solution));
}

template struct ElasticTangentStiffnessModel<2>;
template struct ElasticTangentStiffnessModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
