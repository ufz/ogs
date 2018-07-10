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
    /*material_state_variables*/,  double const T)
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
        std::exp(-_Q / (MaterialLib::PhysicalConstant::IdealGasConstant * T));
    auto const& deviatoric_matrix = Invariants::deviatoric_projection;

    auto const update_residual = [&](ResidualVectorType& r) {

        auto const s_n1 = deviatoric_matrix * solution;
        // ||s_{n+1}|| = sqrt(2.0 * J2)
        double const norm_s_n1 = std::sqrt(2.0 * Invariants::J2(s_n1));
        r = solution - sigma_try +
            2.0 * b * this->_mp.mu(t, x) * std::pow(norm_s_n1, _n - 1) * s_n1;
    };

    auto const update_jacobian = [&](JacobianMatrix& jacobian) {
        auto const s_n1 = deviatoric_matrix * solution;
        double const norm_s_n1 = std::sqrt(2.0 * Invariants::J2(s_n1));
        jacobian =
            KelvinMatrix::Identity() +
            2.0 * b * this->_mp.mu(t, x) * std::pow(norm_s_n1, _n - 1) *
                (deviatoric_matrix +
                 (_n - 1) * std::pow(norm_s_n1, -2) * s_n1 * s_n1.transpose());
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

    KelvinMatrix tangentStiffness = linear_solver.solve(C);

    return {std::make_tuple(solution, createMaterialStateVariables(),
                            tangentStiffness)};
}

template class CreepBGRa<2>;
template class CreepBGRa<3>;

}  // namespace Creep
}  // namespace Solids
}  // namespace MaterialLib
