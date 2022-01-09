/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <cassert>
#include <ostream>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/NewtonRaphson.h"

namespace ProcessLib::RichardsMechanics
{
template <int DisplacementDim>
struct MicroPorosityStateSpace
{
    double phi_m;
    double e_sw;
    double p_L_m;
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> sigma_sw;

    MicroPorosityStateSpace& operator+=(MicroPorosityStateSpace const& state)
    {
        phi_m += state.phi_m;
        e_sw += state.e_sw;
        p_L_m += state.p_L_m;
        sigma_sw += state.sigma_sw;
        return *this;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <int DisplacementDim>
std::ostream& operator<<(std::ostream& os,
                         MicroPorosityStateSpace<DisplacementDim> const& state)
{
    return os << "phi_m: " << state.phi_m << ", "
              << "e_sw: " << state.e_sw << ", "
              << "p_L_m: " << state.p_L_m << ", "
              << "sigma_sw: " << state.sigma_sw.transpose();
}

struct MicroPorosityParameters
{
    NumLib::NewtonRaphsonSolverParameters nonlinear_solver_parameters;

    /// Coefficient \f$\bar\alpha\f$ of the micro structure mass exchange. If
    /// this is given then the micro_saturation MPL property must be provided
    /// too.
    double mass_exchange_coefficient;
};

template <int DisplacementDim>
MicroPorosityStateSpace<DisplacementDim> computeMicroPorosity(
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const&
        I_2_C_el_inverse,
    double const rho_LR_m,  // for simplification equal to rho_LR_M
    double const mu_LR,
    MicroPorosityParameters const& micro_porosity_parameters,
    double const alpha_B, double const phi, double const p_L,
    double const p_L_m_prev,
    MaterialPropertyLib::VariableArray const& /*variables_prev*/,
    double const S_L_m_prev, double const phi_m_prev,
    ParameterLib::SpatialPosition const pos, double const t, double const dt,
    MaterialPropertyLib::Property const& saturation_micro,
    MaterialPropertyLib::Property const& swelling_stress_rate)
{
    namespace MPL = MaterialPropertyLib;
    static constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    // phi_m, e_sw, p_L_m, sigma_sw
    static constexpr int nls_size = 1 + 1 + 1 + kelvin_vector_size;

    static constexpr int i_phi_m = 0;
    static constexpr int i_e_sw = 1;
    static constexpr int i_p_L_m = 2;
    static constexpr int i_sigma_sw = 3;

    using ResidualVectorType = Eigen::Matrix<double, nls_size, 1>;
    using JacobianMatrix =
        Eigen::Matrix<double, nls_size, nls_size, Eigen::RowMajor>;

    Eigen::FullPivLU<Eigen::Matrix<double, nls_size, nls_size, Eigen::RowMajor>>
        linear_solver;

    JacobianMatrix jacobian;

    // Agglomerated solution vector construction.
    ResidualVectorType solution = ResidualVectorType::Zero();

    if (p_L >= 0 && p_L_m_prev >= 0)
    {
        return {solution[i_phi_m], solution[i_e_sw], solution[i_p_L_m],
                solution.template segment<kelvin_vector_size>(i_sigma_sw)};
    }

    double const alpha_bar =
        micro_porosity_parameters.mass_exchange_coefficient;

    auto const update_residual = [&](ResidualVectorType& residual)
    {
        double const delta_phi_m = solution[i_phi_m];
        double const delta_e_sw = solution[i_e_sw];
        auto const& delta_sigma_sw =
            solution.template segment<kelvin_vector_size>(i_sigma_sw);
        double const delta_p_L_m = solution[i_p_L_m];

        double const phi_m = phi_m_prev + delta_phi_m;
        double const p_L_m = p_L_m_prev + delta_p_L_m;
        MPL::VariableArray variables_prev;
        variables_prev[static_cast<int>(MPL::Variable::capillary_pressure)] =
            -p_L_m_prev;
        variables_prev[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L_m_prev;

        MPL::VariableArray variables;
        variables[static_cast<int>(MPL::Variable::capillary_pressure)] = -p_L_m;

        double const S_L_m =
            saturation_micro.template value<double>(variables, pos, t, dt);
        variables[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L_m;
        double const delta_S_L_m = S_L_m - S_L_m_prev;

        auto const sigma_sw_dot =
            MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                swelling_stress_rate.template value<Eigen::Matrix3d>(
                    variables, variables_prev, pos, t, dt));

        residual[i_phi_m] = delta_phi_m - (alpha_B - phi) * delta_e_sw;
        residual[i_e_sw] = delta_e_sw + I_2_C_el_inverse.dot(delta_sigma_sw);
        residual.template segment<kelvin_vector_size>(i_sigma_sw).noalias() =
            delta_sigma_sw - sigma_sw_dot * dt;

        residual[i_p_L_m] =
            rho_LR_m *
                (phi_m * delta_S_L_m - (alpha_B - phi) * S_L_m * delta_e_sw) +
            phi_m * S_L_m * rho_LR_m * delta_e_sw -
            micro_porosity_parameters.mass_exchange_coefficient / mu_LR *
                (p_L - p_L_m) * dt;
    };

    auto const update_jacobian = [&](JacobianMatrix& jacobian)
    {
        jacobian = JacobianMatrix::Identity();

        double const delta_phi_m = solution[i_phi_m];
        double const delta_e_sw = solution[i_e_sw];
        double const delta_p_L_m = solution[i_p_L_m];

        double const phi_m = phi_m_prev + delta_phi_m;
        double const p_L_m = p_L_m_prev + delta_p_L_m;
        MPL::VariableArray variables_prev;
        variables_prev[static_cast<int>(MPL::Variable::capillary_pressure)] =
            -p_L_m_prev;
        MPL::VariableArray variables;
        variables[static_cast<int>(MPL::Variable::capillary_pressure)] = -p_L_m;

        double const S_L_m =
            saturation_micro.template value<double>(variables, pos, t, dt);
        variables_prev[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L_m_prev;
        variables[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L_m;
        double const delta_S_L_m = S_L_m - S_L_m_prev;

        double const dS_L_m_dp_cap_m = saturation_micro.template dValue<double>(
            variables, MPL::Variable::capillary_pressure, pos, t, dt);
        auto const dsigma_sw_dS_L_m =
            MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                swelling_stress_rate.template dValue<Eigen::Matrix3d>(
                    variables, variables_prev, MPL::Variable::liquid_saturation,
                    pos, t, dt));

        jacobian(i_phi_m, i_e_sw) = -(alpha_B - phi);

        jacobian.template block<1, kelvin_vector_size>(i_e_sw, i_sigma_sw) =
            I_2_C_el_inverse.transpose();

        jacobian.template block<kelvin_vector_size, 1>(i_sigma_sw, i_p_L_m) =
            -dsigma_sw_dS_L_m * dS_L_m_dp_cap_m;

        jacobian(i_p_L_m, i_phi_m) =
            rho_LR_m * (delta_S_L_m + S_L_m * delta_e_sw);

        jacobian(i_p_L_m, i_e_sw) = -rho_LR_m * S_L_m * (alpha_B - phi - phi_m);

        jacobian(i_p_L_m, i_p_L_m) =
            alpha_bar / mu_LR * dt -
            rho_LR_m * (phi_m - (alpha_B - phi - phi_m) * delta_e_sw) *
                dS_L_m_dp_cap_m;
    };

    auto const update_solution = [&](ResidualVectorType const& increment)
    { solution += increment; };

    auto newton_solver =
        NumLib::NewtonRaphson<decltype(linear_solver), JacobianMatrix,
                              decltype(update_jacobian), ResidualVectorType,
                              decltype(update_residual),
                              decltype(update_solution)>(
            linear_solver, update_jacobian, update_residual, update_solution,
            micro_porosity_parameters.nonlinear_solver_parameters);

    auto const success_iterations = newton_solver.solve(jacobian);

    if (!success_iterations)
    {
        OGS_FATAL(
            "Could not find solution for local double structure nonlinear "
            "problem.");
    }

    return {solution[i_phi_m], solution[i_e_sw], solution[i_p_L_m],
            solution.template segment<kelvin_vector_size>(i_sigma_sw)};
}
}  // namespace ProcessLib::RichardsMechanics
