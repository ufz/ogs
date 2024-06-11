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

#include <memory>

#include "ConstitutiveRelations/Base.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MathLib/KelvinVector.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveStress_StrainTemperature/SolidMechanics.h"

namespace ProcessLib
{
namespace RichardsMechanics
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatricesTypePressure, int DisplacementDim, int NPoints>
struct IntegrationPointData final
{
    typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatrixTypeDisplacement::GlobalDimNodalMatrixType dNdx_u;

    typename ShapeMatricesTypePressure::NodalRowVectorType N_p;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx_p;

    double integration_weight = std::numeric_limits<double>::quiet_NaN();

    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>
    computeElasticTangentStiffness(
        double const t,
        ParameterLib::SpatialPosition const& x_position,
        double const dt,
        double const temperature,
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material,
        typename MaterialLib::Solids::MechanicsBase<DisplacementDim>::
            MaterialStateVariables const& material_state_variables)
    {
        namespace MPL = MaterialPropertyLib;

        MPL::VariableArray variable_array;
        MPL::VariableArray variable_array_prev;

        using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

        variable_array.stress.emplace<KV>(KV::Zero());
        variable_array.mechanical_strain.emplace<KV>(KV::Zero());
        variable_array.temperature = temperature;

        variable_array_prev.stress.emplace<KV>(KV::Zero());
        variable_array_prev.mechanical_strain.emplace<KV>(KV::Zero());
        variable_array_prev.temperature = temperature;

        auto&& solution = solid_material.integrateStress(
            variable_array_prev, variable_array, t, x_position, dt,
            material_state_variables);

        if (!solution)
        {
            OGS_FATAL("Computation of elastic tangent stiffness failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C =
            std::move(std::get<2>(*solution));

        return C;
    }

    typename BMatricesType::KelvinMatrixType updateConstitutiveRelation(
        MaterialPropertyLib::VariableArray const& variable_array,
        double const t,
        ParameterLib::SpatialPosition const& x_position,
        double const dt,
        double const temperature,
        ProcessLib::ThermoRichardsMechanics::
            ConstitutiveStress_StrainTemperature::EffectiveStressData<
                DisplacementDim>& sigma_eff,
        PrevState<ProcessLib::ThermoRichardsMechanics::
                      ConstitutiveStress_StrainTemperature::EffectiveStressData<
                          DisplacementDim>> const& sigma_eff_prev,
        ProcessLib::ThermoRichardsMechanics::
            ConstitutiveStress_StrainTemperature::MechanicalStrainData<
                DisplacementDim> const&
        /*eps_m*/,
        PrevState<
            ProcessLib::ThermoRichardsMechanics::
                ConstitutiveStress_StrainTemperature::MechanicalStrainData<
                    DisplacementDim>> const& eps_m_prev,
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material,
        std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
            DisplacementDim>::MaterialStateVariables>& material_state_variables)
    {
        MaterialPropertyLib::VariableArray variable_array_prev;
        variable_array_prev.stress = sigma_eff_prev->sigma_eff;
        variable_array_prev.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m_prev->eps_m);
        variable_array_prev.temperature = temperature;

        auto&& solution = solid_material.integrateStress(
            variable_array_prev, variable_array, t, x_position, dt,
            *material_state_variables);

        if (!solution)
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff.sigma_eff, material_state_variables, C) =
            std::move(*solution);

        return C;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace RichardsMechanics
}  // namespace ProcessLib
