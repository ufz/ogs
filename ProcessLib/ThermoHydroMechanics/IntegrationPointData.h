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

#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatricesTypePressure, int DisplacementDim, int NPoints>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
        // Initialize current time step values
        static const int kelvin_vector_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
        sigma_eff.setZero(kelvin_vector_size);
        sigma_eff_ice.setZero(kelvin_vector_size);
        eps.setZero(kelvin_vector_size);
        eps0.setZero(kelvin_vector_size);
        eps0_prev.setZero(kelvin_vector_size);
        eps_m.setZero(kelvin_vector_size);
        eps_m_prev.resize(kelvin_vector_size);
        eps_m_ice.setZero(kelvin_vector_size);
        eps_m_ice_prev.resize(kelvin_vector_size);

        // Previous time step values are not initialized and are set later.
        sigma_eff_prev.resize(kelvin_vector_size);
        sigma_eff_ice_prev.resize(kelvin_vector_size);
    }

    typename BMatricesType::KelvinVectorType sigma_eff, sigma_eff_prev;
    typename BMatricesType::KelvinVectorType eps, eps0, eps0_prev;
    typename BMatricesType::KelvinVectorType eps_m, eps_m_prev;

    typename BMatricesType::KelvinVectorType sigma_eff_ice, sigma_eff_ice_prev;
    typename BMatricesType::KelvinVectorType eps_m_ice, eps_m_ice_prev;

    typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatrixTypeDisplacement::GlobalDimNodalMatrixType dNdx_u;

    // Scalar shape matrices for pressure and temperature.
    typename ShapeMatricesTypePressure::NodalRowVectorType N;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx;

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    double phi_fr = std::numeric_limits<double>::quiet_NaN();
    double phi_fr_prev = std::numeric_limits<double>::quiet_NaN();
    double integration_weight;

    double porosity;

    void pushBackState()
    {
        eps0_prev = eps0;
        eps_m_prev = eps_m;
        sigma_eff_prev = sigma_eff;
        sigma_eff_ice_prev = sigma_eff_ice;
        eps_m_ice_prev = eps_m_ice;
        material_state_variables->pushBackState();
    }

    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>
    computeElasticTangentStiffness(
        double const t,
        ParameterLib::SpatialPosition const& x_position,
        double const dt,
        double const temperature)
    {
        namespace MPL = MaterialPropertyLib;

        MPL::VariableArray variable_array;
        MPL::VariableArray variable_array_prev;

        auto const null_state = solid_material.createMaterialStateVariables();
        solid_material.initializeInternalStateVariables(t, x_position,
                                                        *null_state);

        using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

        variable_array.stress.emplace<KV>(KV::Zero());
        variable_array.mechanical_strain.emplace<KV>(KV::Zero());
        variable_array.temperature = temperature;

        variable_array_prev.stress.emplace<KV>(KV::Zero());
        variable_array_prev.mechanical_strain.emplace<KV>(KV::Zero());
        variable_array_prev.temperature = temperature;

        auto&& solution =
            solid_material.integrateStress(variable_array_prev, variable_array,
                                           t, x_position, dt, *null_state);

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
        double const temperature_prev)
    {
        MaterialPropertyLib::VariableArray variable_array_prev;
        variable_array_prev.stress
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                sigma_eff_prev);
        variable_array_prev.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m_prev);
        variable_array_prev.temperature = temperature_prev;

        auto&& solution = solid_material.integrateStress(
            variable_array_prev, variable_array, t, x_position, dt,
            *material_state_variables);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff, material_state_variables, C) = std::move(*solution);

        return C;
    }

    typename BMatricesType::KelvinMatrixType updateConstitutiveRelationIce(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            ice_constitutive_relation,
        MaterialPropertyLib::VariableArray const& variable_array,
        double const t,
        ParameterLib::SpatialPosition const& x_position,
        double const dt,
        double const temperature_prev)
    {
        MaterialPropertyLib::VariableArray variable_array_prev;

        variable_array_prev.stress
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                sigma_eff_ice_prev);
        variable_array_prev.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m_ice_prev);
        variable_array_prev.temperature = temperature_prev;

        // Extend for non-linear ice materials if necessary.
        auto const null_state =
            ice_constitutive_relation.createMaterialStateVariables();
        ice_constitutive_relation.initializeInternalStateVariables(
            t, x_position, *null_state);
        auto&& solution = ice_constitutive_relation.integrateStress(
            variable_array_prev, variable_array, t, x_position, dt,
            *null_state);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C_IR;
        std::tie(sigma_eff_ice, material_state_variables, C_IR) =
            std::move(*solution);

        return C_IR;
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <int DisplacementDim>
struct IntegrationPointDataForOutput
{
    using DimVector = MatrixPolicyType::VectorType<DisplacementDim>;
    // Darcy velocity for output. Care must be taken for the deactivated
    // elements.
    DimVector velocity = DimVector::Constant(
        DisplacementDim, std::numeric_limits<double>::quiet_NaN());

    double fluid_density = std::numeric_limits<double>::quiet_NaN();
    double viscosity = std::numeric_limits<double>::quiet_NaN();
};

template <int DisplacementDim>
struct ConstitutiveRelationsValues
{
    using DimMatrix =
        typename MatrixPolicyType::MatrixType<DisplacementDim, DisplacementDim>;

    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
        solid_linear_thermal_expansion_coefficient;
    DimMatrix K_over_mu;
    DimMatrix K_pT_thermal_osmosis;
    DimMatrix effective_thermal_conductivity;
    double alpha_biot;
    double beta;
    double beta_SR;
    double c_f;
    double effective_volumetric_heat_capacity;
    double fluid_compressibility;
    double rho;

    // Freezing related values.
    double J_TT_fr;
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> J_uu_fr;
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> J_uT_fr;
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> r_u_fr;
};

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
