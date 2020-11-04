/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace TH2M
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatricesTypePressure, int DisplacementDim, int NPoints>
struct IntegrationPointData final
{
    using GlobalDimMatrixType =
        typename ShapeMatricesTypePressure::GlobalDimMatrixType;
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
        eps.setZero(kelvin_vector_size);
        eps_m.setZero(kelvin_vector_size);
        eps_m_prev.resize(kelvin_vector_size);

        // Previous time step values are not initialized and are set later.
        eps_prev.resize(kelvin_vector_size);
        sigma_eff_prev.resize(kelvin_vector_size);
    }

    typename ShapeMatrixTypeDisplacement::template MatrixType<
        DisplacementDim, NPoints * DisplacementDim>
        N_u_op;
    typename BMatricesType::KelvinVectorType sigma_eff, sigma_eff_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;
    typename BMatricesType::KelvinVectorType eps_m, eps_m_prev;

    typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatrixTypeDisplacement::GlobalDimNodalMatrixType dNdx_u;

    typename ShapeMatricesTypePressure::NodalRowVectorType N_p;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx_p;

    double s_L = std::numeric_limits<double>::quiet_NaN();
    double s_L_prev = std::numeric_limits<double>::quiet_NaN();
    double dsLdpCap = std::numeric_limits<double>::quiet_NaN();
    double MC = std::numeric_limits<double>::quiet_NaN();
    double MW = std::numeric_limits<double>::quiet_NaN();
    double MG = std::numeric_limits<double>::quiet_NaN();

    // phase intrinsic densities
    double rhoGR = std::numeric_limits<double>::quiet_NaN();
    double rhoLR = std::numeric_limits<double>::quiet_NaN();
    double rhoSR = std::numeric_limits<double>::quiet_NaN();
    // solid phase pressure
    double p_SR = std::numeric_limits<double>::quiet_NaN();
    double pWGR = std::numeric_limits<double>::quiet_NaN();

    double dp_vap_dT = std::numeric_limits<double>::quiet_NaN();
    double dp_vap_dpCap = std::numeric_limits<double>::quiet_NaN();

    // real constitutent partial densities
    double rhoCGR = std::numeric_limits<double>::quiet_NaN();
    double rhoCGR_prev = std::numeric_limits<double>::quiet_NaN();
    double rhoWGR = std::numeric_limits<double>::quiet_NaN();
    double rhoWGR_prev = std::numeric_limits<double>::quiet_NaN();
    double rhoCLR = std::numeric_limits<double>::quiet_NaN();
    double rhoCLR_prev = std::numeric_limits<double>::quiet_NaN();
    double rhoWLR = std::numeric_limits<double>::quiet_NaN();
    double rhoWLR_prev = std::numeric_limits<double>::quiet_NaN();

    double drhoGR_dpGR = std::numeric_limits<double>::quiet_NaN();
    double drhoGR_dT = std::numeric_limits<double>::quiet_NaN();

    // phase composition
    // molar fraction
    double xnCG = std::numeric_limits<double>::quiet_NaN();
    double xnWG = std::numeric_limits<double>::quiet_NaN();
    double xnCL = std::numeric_limits<double>::quiet_NaN();
    double xnWL = std::numeric_limits<double>::quiet_NaN();

    double dxnCG_dpGR = std::numeric_limits<double>::quiet_NaN();
    double dxnCG_dpCap = std::numeric_limits<double>::quiet_NaN();
    double dxnCG_dT = std::numeric_limits<double>::quiet_NaN();

    // mass fraction
    double xmCG = std::numeric_limits<double>::quiet_NaN();
    double xmWG = std::numeric_limits<double>::quiet_NaN();
    double xmCL = std::numeric_limits<double>::quiet_NaN();
    double xmWL = std::numeric_limits<double>::quiet_NaN();
    // mass fraction derivatives
    double dxmCG_dpGR = std::numeric_limits<double>::quiet_NaN();
    double dxmWG_dpGR = std::numeric_limits<double>::quiet_NaN();
    double dxmCL_dpLR = std::numeric_limits<double>::quiet_NaN();
    double dxmWL_dpLR = std::numeric_limits<double>::quiet_NaN();
    double dxmCG_dT = std::numeric_limits<double>::quiet_NaN();
    double dxmWG_dT = std::numeric_limits<double>::quiet_NaN();
    double dxmCL_dT = std::numeric_limits<double>::quiet_NaN();
    double dxmWL_dT = std::numeric_limits<double>::quiet_NaN();

    // diffusion coefficients
    double diffusion_coefficient_vapour =
        std::numeric_limits<double>::quiet_NaN();
    double diffusion_coefficient_solvate =
        std::numeric_limits<double>::quiet_NaN();

    // phase enthalpies
    double rho_G_h_G = std::numeric_limits<double>::quiet_NaN();
    double rho_G_h_G_prev = std::numeric_limits<double>::quiet_NaN();
    double rho_L_h_L = std::numeric_limits<double>::quiet_NaN();
    double rho_L_h_L_prev = std::numeric_limits<double>::quiet_NaN();
    double rho_S_h_S = std::numeric_limits<double>::quiet_NaN();
    double rho_S_h_S_prev = std::numeric_limits<double>::quiet_NaN();

    // specific enthalpies
    double h_G = std::numeric_limits<double>::quiet_NaN();
    double h_L = std::numeric_limits<double>::quiet_NaN();
    double h_S = std::numeric_limits<double>::quiet_NaN();
    double h_CG = std::numeric_limits<double>::quiet_NaN();
    double h_WG = std::numeric_limits<double>::quiet_NaN();

    double dhG_dT = std::numeric_limits<double>::quiet_NaN();

    // internal energies
    double u_G = std::numeric_limits<double>::quiet_NaN();
    double u_L = std::numeric_limits<double>::quiet_NaN();
    double rho_u_eff = std::numeric_limits<double>::quiet_NaN();
    double rho_u_eff_prev = std::numeric_limits<double>::quiet_NaN();

    // heat capacities
    double cp_G = std::numeric_limits<double>::quiet_NaN();
    double cp_CG = std::numeric_limits<double>::quiet_NaN();
    double cp_WG = std::numeric_limits<double>::quiet_NaN();
    double cp_L = std::numeric_limits<double>::quiet_NaN();
    double cp_CL = std::numeric_limits<double>::quiet_NaN();
    double cp_WL = std::numeric_limits<double>::quiet_NaN();
    double cp_S = std::numeric_limits<double>::quiet_NaN();

    // porosity
    double phi = std::numeric_limits<double>::quiet_NaN();
    double phi_prev = std::numeric_limits<double>::quiet_NaN();

    double muGR = std::numeric_limits<double>::quiet_NaN();
    double muLR = std::numeric_limits<double>::quiet_NaN();

    GlobalDimMatrixType lambda;

    // other
    double phi_S_p_SR = std::numeric_limits<double>::quiet_NaN();
    double phi_S_p_SR_prev = std::numeric_limits<double>::quiet_NaN();

    double thermal_volume_strain = std::numeric_limits<double>::quiet_NaN();
    double beta_pS = std::numeric_limits<double>::quiet_NaN();
    double beta_T_SR = std::numeric_limits<double>::quiet_NaN();
    double c_p_S = std::numeric_limits<double>::quiet_NaN();
    double alpha_B = std::numeric_limits<double>::quiet_NaN();
    double beta_p_SR = std::numeric_limits<double>::quiet_NaN();

    double k_rel_L = std::numeric_limits<double>::quiet_NaN();
    double k_rel_G = std::numeric_limits<double>::quiet_NaN();

    // solid phase linear thermal expansivity
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> alpha_T_SR;
    GlobalDimMatrixType k_S;

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;
    double integration_weight = std::numeric_limits<double>::quiet_NaN();

    void pushBackState()
    {
        eps_m_prev = eps_m;
        sigma_eff_prev = sigma_eff;
        s_L_prev = s_L;

        rho_G_h_G_prev = rho_G_h_G;
        rho_L_h_L_prev = rho_L_h_L;
        rho_S_h_S_prev = rho_S_h_S;

        phi_prev = phi;

        phi_S_p_SR_prev = phi_S_p_SR;

        rhoCGR_prev = rhoCGR;
        rhoWGR_prev = rhoWGR;
        rhoCLR_prev = rhoCLR;
        rhoWLR_prev = rhoWLR;

        rho_u_eff_prev = rho_u_eff;

        material_state_variables->pushBackState();
    }

    template <typename DisplacementVectorType>
    typename BMatricesType::KelvinMatrixType updateConstitutiveRelation(
        double const t,
        ParameterLib::SpatialPosition const& x_position,
        double const dt,
        DisplacementVectorType const& /*u*/,
        double const T)
    {
        // TODO (Norbert) These current time step variables should be filled in
        // the FEM and passed here instead of the T and thermal_strain.
        MaterialPropertyLib::VariableArray variable_array;
        variable_array[static_cast<int>(MaterialPropertyLib::Variable::stress)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                sigma_eff);
        variable_array[static_cast<int>(
                           MaterialPropertyLib::Variable::mechanical_strain)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);
        variable_array[static_cast<int>(
                           MaterialPropertyLib::Variable::temperature)]
            .emplace<double>(T);

        MaterialPropertyLib::VariableArray variable_array_prev;
        variable_array_prev[static_cast<int>(
                                MaterialPropertyLib::Variable::stress)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                sigma_eff_prev);
        variable_array_prev[static_cast<int>(MaterialPropertyLib::Variable::
                                                 mechanical_strain)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m_prev);
        variable_array_prev[static_cast<int>(
                                MaterialPropertyLib::Variable::temperature)]
            .emplace<double>(T);  // TODO (Norbert) should be T_prev
        auto&& solution = solid_material.integrateStress(
            variable_array_prev, variable_array, t, x_position, dt,
            *material_state_variables);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff, material_state_variables, C) = std::move(*solution);

        return C;
    }

    template <typename DisplacementVectorType>
    typename BMatricesType::KelvinMatrixType updateConstitutiveRelationThermal(
        double const t,
        ParameterLib::SpatialPosition const& x_position,
        double const dt,
        DisplacementVectorType const& /*u*/,
        double const T,
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const&
            thermal_strain)
    {
        eps_m.noalias() = eps - thermal_strain;

        // TODO (Norbert) These current time step variables should be filled in
        // the FEM and passed here instead of the T and thermal_strain.
        MaterialPropertyLib::VariableArray variable_array;
        variable_array[static_cast<int>(MaterialPropertyLib::Variable::stress)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                sigma_eff);
        variable_array[static_cast<int>(
                           MaterialPropertyLib::Variable::mechanical_strain)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);
        variable_array[static_cast<int>(
                           MaterialPropertyLib::Variable::temperature)]
            .emplace<double>(T);

        MaterialPropertyLib::VariableArray variable_array_prev;
        variable_array_prev[static_cast<int>(
                                MaterialPropertyLib::Variable::stress)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                sigma_eff_prev);
        variable_array_prev[static_cast<int>(MaterialPropertyLib::Variable::
                                                 mechanical_strain)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m_prev);
        variable_array_prev[static_cast<int>(
                                MaterialPropertyLib::Variable::temperature)]
            .emplace<double>(T);  // TODO (Norbert) should be T_prev

        auto&& solution = solid_material.integrateStress(
            variable_array_prev, variable_array, t, x_position, dt,
            *material_state_variables);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff, material_state_variables, C) = std::move(*solution);

        return C;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace TH2M
}  // namespace ProcessLib
