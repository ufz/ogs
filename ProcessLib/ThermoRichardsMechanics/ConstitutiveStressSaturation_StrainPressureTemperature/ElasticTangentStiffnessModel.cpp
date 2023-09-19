/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElasticTangentStiffnessModel.h"

#include "MaterialLib/SolidModels/MFront/Variable.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
template <int DisplacementDim>
void ElasticTangentStiffnessModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, TemperatureData<DisplacementDim> const& T_data,
    ElasticTangentStiffnessData<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;

    auto const null_state = solid_material_.createMaterialStateVariables();
    solid_material_.initializeInternalStateVariables(x_t.t, x_t.x, *null_state);

    using KV = KelvinVector<DisplacementDim>;

    MPL::VariableArray variable_array;
    {
        // gradients
        variable_array.mechanical_strain = KV::Zero().eval();
        variable_array.liquid_phase_pressure = 0.0;

        // external state variables
        variable_array.temperature = T_data.T;
    }

    MPL::VariableArray variable_array_prev;
    {
        // thermodynamic forces
        variable_array_prev.stress = KV::Zero().eval();
        variable_array_prev.liquid_saturation = 1.0;

        // gradients
        variable_array_prev.mechanical_strain = KV::Zero().eval();
        variable_array_prev.liquid_phase_pressure = 0.0;

        // external state variables
        variable_array_prev.temperature = T_data.T_prev;
    }

    auto&& solution = solid_material_.integrateStress(
        variable_array_prev, variable_array, x_t.t, x_t.x, x_t.dt, *null_state);

    if (!solution)
    {
        OGS_FATAL("Computation of elastic tangent stiffness failed.");
    }

    auto const& tangent_operator_data = std::get<2>(*solution);
    out.C_el = tangent_operator_blocks_view_.block(MSM::stress, MSM::strain,
                                                   tangent_operator_data);
}

template struct ElasticTangentStiffnessModel<2>;
template struct ElasticTangentStiffnessModel<3>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
