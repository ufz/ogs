// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <limits>
#include <memory>

#include "MaterialLib/MPL/VariableType.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{
using SpaceTimeData = ConstitutiveRelations::SpaceTimeData;

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
        // Initialise current time step values
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

        // Previous time step values are not initialised and are set later.
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
    double integration_weight = 0;

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
    computeElasticTangentStiffness(SpaceTimeData const& space_time_data,
                                   double const temperature)
    {
        auto const& [x_position, t, dt] = space_time_data;
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
            *material_state_variables);

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
        SpaceTimeData const& space_time_data, double const temperature_prev)
    {
        auto const& [x_position, t, dt] = space_time_data;
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

    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>
    computeElasticTangentStiffnessIce(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            ice_constitutive_relation,
        SpaceTimeData const& space_time_data, double const temperature)
    {
        auto const& [x_position, t, dt] = space_time_data;
        namespace MPL = MaterialPropertyLib;

        MPL::VariableArray variable_array;
        MPL::VariableArray variable_array_prev;

        auto const null_state =
            ice_constitutive_relation.createMaterialStateVariables();
        ice_constitutive_relation.initializeInternalStateVariables(
            t, x_position, *null_state);

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

    typename BMatricesType::KelvinMatrixType updateConstitutiveRelationIce(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            ice_constitutive_relation,
        MaterialPropertyLib::VariableArray const& variable_array,
        SpaceTimeData const& space_time_data, double const temperature_prev)
    {
        auto const& [x_position, t, dt] = space_time_data;
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
    double rho_fr = std::numeric_limits<double>::quiet_NaN();
};

/// The integration point state of one property evaluation: the integration
/// point data and the cached output quantities of a single integration point,
/// referenced rather than owned. Evaluated property values are written back
/// to the referenced state as documented for each update*() function.
template <typename IpData, int DisplacementDim>
struct IntegrationPointState final
{
    IpData& ip_data;
    IntegrationPointDataForOutput<DisplacementDim>& ip_data_output;
};

/// Values shared by all balance equations. Evaluated once per integration
/// point by updateCommonProperties(); every other update*Properties() function
/// takes these as input and must not recompute them.
template <int DisplacementDim>
struct CommonPropertyValues final
{
    MaterialPropertyLib::VariableArray vars;

    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> eps_prev =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>::Constant(
            std::numeric_limits<double>::quiet_NaN());
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
        solid_linear_thermal_expansion_coefficient =
            MathLib::KelvinVector::KelvinVectorType<DisplacementDim>::Constant(
                std::numeric_limits<double>::quiet_NaN());
    double solid_volumetric_thermal_expansion_coefficient =
        std::numeric_limits<double>::quiet_NaN();

    double T_prev_int_pt = std::numeric_limits<double>::quiet_NaN();
    double dT_int_pt = std::numeric_limits<double>::quiet_NaN();
    double T_dot = std::numeric_limits<double>::quiet_NaN();
    double p_prev_int_pt = std::numeric_limits<double>::quiet_NaN();
    double p_dot = std::numeric_limits<double>::quiet_NaN();
    double eps_v_dot = std::numeric_limits<double>::quiet_NaN();

    double alpha_biot = std::numeric_limits<double>::quiet_NaN();
    /// Biot coefficient including the ice volume fraction correction. Equal to
    /// alpha_biot without a frozen liquid phase.
    double modified_biot = std::numeric_limits<double>::quiet_NaN();
    double drho_LR_dT = std::numeric_limits<double>::quiet_NaN();
    double solid_density = std::numeric_limits<double>::quiet_NaN();

    /// \name Freezing related values; NaN without a frozen liquid phase.
    /// @{
    double dphi_fr_dT = std::numeric_limits<double>::quiet_NaN();
    double d2phi_fr_dT2 = std::numeric_limits<double>::quiet_NaN();
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
        ice_linear_thermal_expansion_coefficient =
            MathLib::KelvinVector::KelvinVectorType<DisplacementDim>::Constant(
                std::numeric_limits<double>::quiet_NaN());
    /// @}
};

/// Darcy flow related values needed by the energy and mass balance equations.
template <int DisplacementDim>
struct FlowPropertyValues final
{
    using DimVector = MatrixPolicyType::VectorType<DisplacementDim>;
    using DimMatrix =
        typename MatrixPolicyType::MatrixType<DisplacementDim, DisplacementDim>;

    DimMatrix K_over_mu =
        DimMatrix::Constant(DisplacementDim, DisplacementDim,
                            std::numeric_limits<double>::quiet_NaN());
    DimMatrix K_pT_thermal_osmosis =
        DimMatrix::Constant(DisplacementDim, DisplacementDim,
                            std::numeric_limits<double>::quiet_NaN());
    DimVector grad_p = DimVector::Constant(
        DisplacementDim, std::numeric_limits<double>::quiet_NaN());
    double k_rel = std::numeric_limits<double>::quiet_NaN();
    DimVector buoyancy_flux = DimVector::Constant(
        DisplacementDim, std::numeric_limits<double>::quiet_NaN());
};

/// Specific heat capacities and the sensible volumetric heat capacity.
struct HeatCapacityValues final
{
    double c_f = std::numeric_limits<double>::quiet_NaN();
    double c_s = std::numeric_limits<double>::quiet_NaN();
    /// NaN without a frozen liquid phase.
    double c_fr = std::numeric_limits<double>::quiet_NaN();
    /// Without the latent contribution.
    double sensible_volumetric_heat_capacity =
        std::numeric_limits<double>::quiet_NaN();
};

template <int DisplacementDim>
struct ThermalPropertyValues final
{
    using DimMatrix =
        typename MatrixPolicyType::MatrixType<DisplacementDim, DisplacementDim>;
    using DimVector = MatrixPolicyType::VectorType<DisplacementDim>;

    DimMatrix effective_thermal_conductivity =
        DimMatrix::Constant(DisplacementDim, DisplacementDim,
                            std::numeric_limits<double>::quiet_NaN());
    DimMatrix dlambda_eff_dT =
        DimMatrix::Constant(DisplacementDim, DisplacementDim,
                            std::numeric_limits<double>::quiet_NaN());
    DimVector dvelocity_dT = DimVector::Constant(
        DisplacementDim, std::numeric_limits<double>::quiet_NaN());

    double dk_rel_dT = std::numeric_limits<double>::quiet_NaN();
    double c_f = std::numeric_limits<double>::quiet_NaN();
    // Sensible contribution to the volumetric heat capacity. The effective
    // value used in the temperature mass matrix is sensible - latent.
    double sensible_volumetric_heat_capacity =
        std::numeric_limits<double>::quiet_NaN();
    // Latent contribution L (effective = sensible - L). Zero without a
    // frozen liquid phase.
    double latent_volumetric_heat_capacity = 0.0;

    double J_TT = std::numeric_limits<double>::quiet_NaN();
};

struct HydraulicPropertyValues final
{
    double drho_L_dp = std::numeric_limits<double>::quiet_NaN();
    /* // TODO(WW): Add this second derivative in the future once the
      computation of the secondary derivatives of the corresponding MPL
      properties are available.
        double d2rho_LR_dp2;
     */
    double solid_skeleton_compressibility =
        std::numeric_limits<double>::quiet_NaN();
    double fluid_compressibility = std::numeric_limits<double>::quiet_NaN();
    double beta = std::numeric_limits<double>::quiet_NaN();
    double beta_SR = std::numeric_limits<double>::quiet_NaN();

    /// \name Freezing related values; NaN without a frozen liquid phase.
    /// @{
    double storage_p_fr = std::numeric_limits<double>::quiet_NaN();
    double storage_T_fr = std::numeric_limits<double>::quiet_NaN();
    double J_pT_fr = std::numeric_limits<double>::quiet_NaN();
    /// @}
};

template <int DisplacementDim>
struct MechanicalPropertyValues final
{
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>::Constant(
            std::numeric_limits<double>::quiet_NaN());
    double rho = std::numeric_limits<double>::quiet_NaN();

    /// \name Freezing related values; unset without a frozen liquid phase.
    /// @{
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> J_uu_fr =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>::Constant(
            std::numeric_limits<double>::quiet_NaN());
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> r_u_fr =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>::Constant(
            std::numeric_limits<double>::quiet_NaN());
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> J_uT_fr =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>::Constant(
            std::numeric_limits<double>::quiet_NaN());
    /// @}
};

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
