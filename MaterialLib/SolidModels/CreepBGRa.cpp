/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreepBGRa.cpp
 *  Created on July 6, 2018, 9:53 AM
 */

#include "CreepBGRa.h"

#include "BaseLib/Error.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialLib
{
namespace Solids
{
namespace Creep
{
template <int DisplacementDim>
boost::optional<std::tuple<typename CreepBGRa<DisplacementDim>::KelvinVector,
                           std::unique_ptr<typename MechanicsBase<
                               DisplacementDim>::MaterialStateVariables>,
                           typename CreepBGRa<DisplacementDim>::KelvinMatrix>>
CreepBGRa<DisplacementDim>::integrateStress(
    double const t, ProcessLib::SpatialPosition const& x, double const dt,
    KelvinVector const& eps_prev, KelvinVector const& eps,
    KelvinVector const& sigma_prev,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
    /*material_state_variables*/,
    double const T)
{
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    Eigen::FullPivLU<Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize,
                                   Eigen::RowMajor>>
        linear_solver;

    const auto C = this->getElasticTensor(t, x, T);
    KelvinVector sigma_try = sigma_prev + C * (eps - eps_prev);
    ResidualVectorType solution = sigma_try;

    const double b =
        dt * _coef *
        std::exp(-_q / (MaterialLib::PhysicalConstant::IdealGasConstant * T));
    auto const& deviatoric_matrix = Invariants::deviatoric_projection;

    // In newton_solver.solve(), the Jacobian is calculated first, and then
    // then comes the assembly of the residue vector. In order to save
    // computation time, the following two variables are set visible in this
    // function. The two variables keep the computation results of the
    // norm of the deviatoric stress and 2bG||s_{n+1}||^{n-1} in
    // update_jacobian, and then they are directly used in update_residual
    // without repeating the computation. Of course, update_jacobian is not a
    // pure function anymore, but has side effects.
    double pow_norm_s_n1_n_minus_one_2b_G = 0.;

    KelvinVector s_n1;
    auto const update_jacobian = [&](JacobianMatrix& jacobian) {
        // side effect
        s_n1 = deviatoric_matrix * solution;
        double const norm_s_n1 = Invariants::FrobeniusNorm(s_n1);
        // side effect
        pow_norm_s_n1_n_minus_one_2b_G =
            2.0 * b * this->_mp.mu(t, x) * std::pow(norm_s_n1, _n - 1);
        jacobian = KelvinMatrix::Identity() +
                   pow_norm_s_n1_n_minus_one_2b_G *
                       (deviatoric_matrix +
                        ((_n - 1) / (norm_s_n1 * norm_s_n1)) * s_n1 *
                            s_n1.transpose());
    };

    auto const update_residual = [&](ResidualVectorType& r) {
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
        return {};

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
    double const t, double const dt, ProcessLib::SpatialPosition const& x,
    double const T, double const deviatoric_stress_norm) const
{
    return 2.0 * _coef *
           std::exp(-_q /
                    (MaterialLib::PhysicalConstant::IdealGasConstant * T)) *
           this->_mp.mu(t, x) * std::pow(deviatoric_stress_norm, _n - 1) * dt *
           _q / (MaterialLib::PhysicalConstant::IdealGasConstant * T * T);
}

template class CreepBGRa<2>;
template class CreepBGRa<3>;

}  // namespace Creep
}  // namespace Solids
}  // namespace MaterialLib
