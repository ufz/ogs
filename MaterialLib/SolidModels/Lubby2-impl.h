/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/NewtonRaphson.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
bool Lubby2<DisplacementDim>::computeConstitutiveRelation(
    double const t,
    ProcessLib::SpatialPosition const& x,
    double const dt,
    KelvinVector const& /*eps_prev*/,
    KelvinVector const& eps,
    KelvinVector const& /*sigma_prev*/,
    KelvinVector& sigma,
    KelvinMatrix& C,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
        material_state_variables)
{
    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;

    assert(dynamic_cast<MaterialStateVariables*>(&material_state_variables) !=
           nullptr);
    MaterialStateVariables& state =
        static_cast<MaterialStateVariables&>(material_state_variables);
    // calculation of deviatoric parts
    auto const& P_dev = Invariants::deviatoric_projection;
    KelvinVector const epsd_i = P_dev * eps;

    // initial guess as elastic predictor
    KelvinVector sigd_j = 2.0 * (epsd_i - state.eps_M_t - state.eps_K_t);

    // Calculate effective stress and update material properties
    double sig_eff = Invariants::equivalentStress(sigd_j);
    updateBurgersProperties(t, x, sig_eff * state.GM, state);

    using LocalJacobianMatrix =
        Eigen::Matrix<double, KelvinVectorSize * 3, KelvinVectorSize * 3,
                      Eigen::RowMajor>;

    // Linear solver for the newton loop is required after the loop with the
    // same matrix. This saves one decomposition.
    Eigen::PartialPivLU<LocalJacobianMatrix> linear_solver(KelvinVectorSize *
                                                           3);

    // Different solvers are available for the solution of the local system.
    // TODO Make the following choice of linear solvers available from the
    // input file configuration:
    //      K_loc.partialPivLu().solve(-res_loc);
    //      K_loc.fullPivLu().solve(-res_loc);
    //      K_loc.householderQr().solve(-res_loc);
    //      K_loc.colPivHouseholderQr().solve(res_loc);
    //      K_loc.fullPivHouseholderQr().solve(-res_loc);
    //      K_loc.llt().solve(-res_loc);
    //      K_loc.ldlt().solve(-res_loc);

    {  // Local Newton solver
        using LocalResidualVector =
            Eigen::Matrix<double, KelvinVectorSize * 3, 1>;

        LocalJacobianMatrix K_loc;
        auto const update_residual = [&](LocalResidualVector& residual) {
            calculateResidualBurgers(dt, epsd_i, sigd_j, state.eps_K_j,
                                     state.eps_K_t, state.eps_M_j,
                                     state.eps_M_t, residual, state);
        };

        auto const update_jacobian = [&](LocalJacobianMatrix& jacobian) {
            calculateJacobianBurgers(
                t, x, dt, jacobian, sig_eff, sigd_j, state.eps_K_j,
                state);  // for solution dependent Jacobians
        };

        auto const update_solution = [&](LocalResidualVector const& increment) {
            // increment solution vectors
            sigd_j.noalias() += increment.template block<KelvinVectorSize, 1>(
                KelvinVectorSize * 0, 0);
            state.eps_K_j.noalias() +=
                increment.template block<KelvinVectorSize, 1>(
                    KelvinVectorSize * 1, 0);
            state.eps_M_j.noalias() +=
                increment.template block<KelvinVectorSize, 1>(
                    KelvinVectorSize * 2, 0);

            // Calculate effective stress and update material properties
            sig_eff = MaterialLib::SolidModels::Invariants<
                KelvinVectorSize>::equivalentStress(sigd_j);
            updateBurgersProperties(t, x, sig_eff * state.GM, state);
        };

        // TODO Make the following choice of maximum iterations and convergence
        // criteria available from the input file configuration:
        const int maximum_iterations(20);
        const double tolerance(1.e-10);

        auto newton_solver = NumLib::NewtonRaphson<
            decltype(linear_solver), LocalJacobianMatrix,
            decltype(update_jacobian), LocalResidualVector,
            decltype(update_residual), decltype(update_solution)>(
            linear_solver, update_jacobian, update_residual, update_solution,
            maximum_iterations, tolerance);

        auto const success_iterations = newton_solver.solve(K_loc);

        if (!success_iterations)
            return false;

        // If the Newton loop didn't run, the linear solver will not be
        // initialized.
        // This happens usually for the first iteration of the first timestep.
        if (*success_iterations == 0)
            linear_solver.compute(K_loc);
    }

    // Hydrostatic part for the stress and the tangent.
    double const eps_i_trace = Invariants::trace(eps);

    sigma.noalias() =
        state.GM * sigd_j + state.KM * eps_i_trace * Invariants::identity2;

    // Calculate dGdE for time step
    Eigen::Matrix<double, KelvinVectorSize * 3, KelvinVectorSize,
                  Eigen::RowMajor> const dGdE = calculatedGdEBurgers();

    // Consistent tangent from local Newton iteration of material
    // functionals.
    // Only the upper left block is relevant for the global tangent.
    auto dzdE = linear_solver.solve(-dGdE)
                    .template block<KelvinVectorSize, KelvinVectorSize>(0, 0);

    auto const& P_sph = Invariants::spherical_projection;
    C.noalias() = state.GM * dzdE * P_dev + 3. * state.KM * P_sph;

    return true;
}

template <int DisplacementDim>
void Lubby2<DisplacementDim>::updateBurgersProperties(
    double const t,
    ProcessLib::SpatialPosition const& x,
    double s_eff,
    MaterialStateVariables& state)
{
    state.GM = _mp.GM0(t, x)[0];
    state.KM = _mp.KM0(t, x)[0];
    state.GK = _mp.GK0(t, x)[0] * std::exp(_mp.mK(t, x)[0] * s_eff);
    state.etaK = _mp.etaK0(t, x)[0] * std::exp(_mp.mvK(t, x)[0] * s_eff);
    state.etaM = _mp.etaM0(t, x)[0] * std::exp(_mp.mvM(t, x)[0] * s_eff);
}

template <int DisplacementDim>
void Lubby2<DisplacementDim>::calculateResidualBurgers(
    const double dt,
    const KelvinVector& strain_curr,
    const KelvinVector& stress_curr,
    KelvinVector& strain_Kel_curr,
    const KelvinVector& strain_Kel_t,
    KelvinVector& strain_Max_curr,
    const KelvinVector& strain_Max_t,
    ResidualVector& res,
    MaterialStateVariables const& state)
{
    // calculate stress residual
    res.template block<KelvinVectorSize, 1>(0, 0).noalias() =
        stress_curr - 2. * (strain_curr - strain_Kel_curr - strain_Max_curr);

    // calculate Kelvin strain residual
    res.template block<KelvinVectorSize, 1>(KelvinVectorSize, 0).noalias() =
        1. / dt * (strain_Kel_curr - strain_Kel_t) -
        1. / (2. * state.etaK) *
            (state.GM * stress_curr - 2. * state.GK * strain_Kel_curr);

    // calculate Maxwell strain residual
    res.template block<KelvinVectorSize, 1>(2 * KelvinVectorSize, 0).noalias() =
        1. / dt * (strain_Max_curr - strain_Max_t) -
        0.5 * state.GM / state.etaM * stress_curr;
}

template <int DisplacementDim>
void Lubby2<DisplacementDim>::calculateJacobianBurgers(
    double const t,
    ProcessLib::SpatialPosition const& x,
    const double dt,
    JacobianMatrix& Jac,
    double s_eff,
    const KelvinVector& sig_i,
    const KelvinVector& eps_K_i,
    MaterialStateVariables const& state)
{
    Jac.setZero();

    // build G_11
    Jac.template block<KelvinVectorSize, KelvinVectorSize>(0, 0).setIdentity();

    // build G_12
    Jac.template block<KelvinVectorSize, KelvinVectorSize>(0, KelvinVectorSize)
        .diagonal()
        .setConstant(2);

    // build G_13
    Jac.template block<KelvinVectorSize, KelvinVectorSize>(0,
                                                           2 * KelvinVectorSize)
        .diagonal()
        .setConstant(2);

    // build G_21
    Jac.template block<KelvinVectorSize, KelvinVectorSize>(KelvinVectorSize, 0)
        .noalias() = -0.5 * state.GM / state.etaK * KelvinMatrix::Identity();
    if (s_eff > 0.)
    {
        KelvinVector const eps_K_aid =
            1. / (state.etaK * state.etaK) *
            (state.GM * sig_i - 2. * state.GK * eps_K_i);

        KelvinVector const dG_K =
            1.5 * _mp.mK(t, x)[0] * state.GK * state.GM / s_eff * sig_i;
        KelvinVector const dmu_vK =
            1.5 * _mp.mvK(t, x)[0] * state.GM * state.etaK / s_eff * sig_i;
        Jac.template block<KelvinVectorSize, KelvinVectorSize>(KelvinVectorSize,
                                                               0)
            .noalias() += 0.5 * eps_K_aid * dmu_vK.transpose() +
                          1. / state.etaK * eps_K_i * dG_K.transpose();
    }

    // build G_22
    Jac.template block<KelvinVectorSize, KelvinVectorSize>(KelvinVectorSize,
                                                           KelvinVectorSize)
        .diagonal()
        .setConstant(1. / dt + state.GK / state.etaK);

    // nothing to do for G_23

    // build G_31
    Jac.template block<KelvinVectorSize, KelvinVectorSize>(2 * KelvinVectorSize,
                                                           0)
        .noalias() = -0.5 * state.GM / state.etaM * KelvinMatrix::Identity();
    if (s_eff > 0.)
    {
        KelvinVector const dmu_vM =
            1.5 * _mp.mvM(t, x)[0] * state.GM * state.etaM / s_eff * sig_i;
        Jac.template block<KelvinVectorSize, KelvinVectorSize>(
               2 * KelvinVectorSize, 0)
            .noalias() += 0.5 * state.GM / (state.etaM * state.etaM) * sig_i *
                          dmu_vM.transpose();
    }

    // nothing to do for G_32

    // build G_33
    Jac.template block<KelvinVectorSize, KelvinVectorSize>(2 * KelvinVectorSize,
                                                           2 * KelvinVectorSize)
        .diagonal()
        .setConstant(1. / dt);
}

}  // namespace Solids
}  // namespace MaterialLib
