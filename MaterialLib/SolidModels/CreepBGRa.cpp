/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on July 6, 2018, 9:53 AM
 */

#include "CreepBGRa.h"

#include <limits>

#include "MaterialLib/MPL/Utils/GetSymmetricTensor.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MPL = MaterialPropertyLib;

namespace MaterialLib
{
namespace Solids
{
namespace Creep
{
double getCreepConstantCoefficient(const double A, const double n,
                                   const double sigma0)
{
    return A * std::pow(1.5, 0.5 * (1 + n)) / std::pow(sigma0, n);
}

template <int DisplacementDim>
std::optional<std::tuple<typename CreepBGRa<DisplacementDim>::KelvinVector,
                         std::unique_ptr<typename MechanicsBase<
                             DisplacementDim>::MaterialStateVariables>,
                         typename CreepBGRa<DisplacementDim>::KelvinMatrix>>
CreepBGRa<DisplacementDim>::integrateStress(
    MaterialPropertyLib::VariableArray const& variable_array_prev,
    MaterialPropertyLib::VariableArray const& variable_array, double const t,
    ParameterLib::SpatialPosition const& x, double const dt,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
    /*material_state_variables*/) const
{
    auto const& eps_m = std::get<MPL::SymmetricTensor<DisplacementDim>>(
        variable_array[static_cast<int>(MPL::Variable::mechanical_strain)]);
    auto const& eps_m_prev = std::get<MPL::SymmetricTensor<DisplacementDim>>(
        variable_array_prev[static_cast<int>(
            MPL::Variable::mechanical_strain)]);
    auto const& sigma_prev = std::get<MPL::SymmetricTensor<DisplacementDim>>(
        variable_array_prev[static_cast<int>(MPL::Variable::stress)]);
    auto const T = std::get<double>(
        variable_array_prev[static_cast<int>(MPL::Variable::temperature)]);

    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    Eigen::FullPivLU<Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize,
                                   Eigen::RowMajor>>
        linear_solver;

    const auto C = this->getElasticTensor(t, x, T);
    KelvinVector sigma_try = sigma_prev + C * (eps_m - eps_m_prev);

    auto const& deviatoric_matrix = Invariants::deviatoric_projection;

    double const norm_s_try =
        Invariants::FrobeniusNorm(deviatoric_matrix * sigma_try);
    // In case |s_{try}| is zero and _n < 3 (rare case).
    if (norm_s_try < std::numeric_limits<double>::epsilon() * C(0, 0))
    {
        return {std::make_tuple(sigma_try, createMaterialStateVariables(), C)};
    }

    ResidualVectorType solution = sigma_try;

    const double A = _a(t, x)[0];
    const double n = _n(t, x)[0];
    const double sigma0 = _sigma_f(t, x)[0];
    const double Q = _q(t, x)[0];

    const double constant_coefficient =
        getCreepConstantCoefficient(A, n, sigma0);

    const double b =
        dt * constant_coefficient *
        std::exp(-Q / (MaterialLib::PhysicalConstant::IdealGasConstant * T));

    double const G2b = 2.0 * b * this->_mp.mu(t, x);

    auto const update_jacobian = [&solution, &G2b,
                                  &n](JacobianMatrix& jacobian) {
        auto const& D = Invariants::deviatoric_projection;
        KelvinVector const s_n1 = D * solution;
        double const norm_s_n1 = Invariants::FrobeniusNorm(s_n1);
        double const pow_norm_s_n1_n_minus_one_2b_G =
            G2b * std::pow(norm_s_n1, n - 1);
        jacobian = KelvinMatrix::Identity() +
                   (pow_norm_s_n1_n_minus_one_2b_G * D +
                    (n - 1) * G2b * std::pow(norm_s_n1, n - 3) * s_n1 *
                        s_n1.transpose());
    };

    auto const update_residual = [&solution, &G2b, &n,
                                  &sigma_try](ResidualVectorType& r) {
        KelvinVector const s_n1 = Invariants::deviatoric_projection * solution;
        double const norm_s_n1 = Invariants::FrobeniusNorm(s_n1);
        double const pow_norm_s_n1_n_minus_one_2b_G =
            G2b * std::pow(norm_s_n1, n - 1);
        r = solution - sigma_try + pow_norm_s_n1_n_minus_one_2b_G * s_n1;
    };

    auto const update_solution = [&](ResidualVectorType const& increment) {
        solution += increment;
    };

    auto newton_solver =
        NumLib::NewtonRaphson<decltype(linear_solver), JacobianMatrix,
                              decltype(update_jacobian), ResidualVectorType,
                              decltype(update_residual),
                              decltype(update_solution)>(
            linear_solver, update_jacobian, update_residual, update_solution,
            _nonlinear_solver_parameters);

    JacobianMatrix jacobian;
    auto const success_iterations = newton_solver.solve(jacobian);

    if (!success_iterations)
    {
        return {};
    }

    // If *success_iterations>0, tangentStiffness = J_(sigma)^{-1}C
    // where J_(sigma) is the Jacobian of the last local Newton-Raphson
    // iteration, which is already LU decomposed.
    KelvinMatrix tangentStiffness =
        (*success_iterations == 0) ? C : linear_solver.solve(C);

    return {std::make_tuple(solution, createMaterialStateVariables(),
                            tangentStiffness)};
}

template <int DisplacementDim>
double CreepBGRa<DisplacementDim>::getTemperatureRelatedCoefficient(
    double const t, double const dt, ParameterLib::SpatialPosition const& x,
    double const T, double const deviatoric_stress_norm) const
{
    const double A = _a(t, x)[0];
    const double n = _n(t, x)[0];
    const double sigma0 = _sigma_f(t, x)[0];
    const double Q = _q(t, x)[0];

    const double constant_coefficient =
        getCreepConstantCoefficient(A, n, sigma0);

    return 2.0 * constant_coefficient *
           std::exp(-Q /
                    (MaterialLib::PhysicalConstant::IdealGasConstant * T)) *
           this->_mp.mu(t, x) * std::pow(deviatoric_stress_norm, n - 1) * dt *
           Q / (MaterialLib::PhysicalConstant::IdealGasConstant * T * T);
}

template class CreepBGRa<2>;
template class CreepBGRa<3>;

}  // namespace Creep
}  // namespace Solids
}  // namespace MaterialLib
