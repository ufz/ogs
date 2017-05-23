/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace MaterialLib
{
namespace Solids
{
namespace Lubby2
{
/// Calculates the 18x6 derivative of the residuals with respect to total
/// strain.
///
/// Function definition can not be moved into implementation because of a
/// MSVC compiler errors. See
/// http://stackoverflow.com/questions/1484885/strange-vc-compile-error-c2244
/// and https://support.microsoft.com/en-us/kb/930198
template <int DisplacementDim>
Eigen::Matrix<double, Lubby2<DisplacementDim>::JacobianResidualSize,
              Lubby2<DisplacementDim>::KelvinVectorSize>
calculatedGdEBurgers()
{
    Eigen::Matrix<double, Lubby2<DisplacementDim>::JacobianResidualSize,
                  Lubby2<DisplacementDim>::KelvinVectorSize>
        dGdE;
    dGdE.setZero();
    dGdE.template block<Lubby2<DisplacementDim>::KelvinVectorSize,
                        Lubby2<DisplacementDim>::KelvinVectorSize>(0, 0)
        .diagonal()
        .setConstant(-2);
    return dGdE;
}

template <int DisplacementDim, typename LinearSolver>
ProcessLib::KelvinMatrixType<DisplacementDim> tangentStiffnessA(
    double const GM0, double const KM0, LinearSolver const& linear_solver)
{
    // Calculate dGdE for time step
    auto const dGdE = calculatedGdEBurgers<DisplacementDim>();

    // Consistent tangent from local Newton iteration of material
    // functionals.
    // Only the upper left block is relevant for the global tangent.
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

    KelvinMatrix const dzdE =
        linear_solver.solve(-dGdE)
            .template block<KelvinVectorSize, KelvinVectorSize>(0, 0);

    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;
    auto const& P_sph = Invariants::spherical_projection;
    auto const& P_dev = Invariants::deviatoric_projection;

    KelvinMatrix C = GM0 * dzdE * P_dev + 3. * KM0 * P_sph;
    return C;
};

template <int DisplacementDim>
boost::optional<std::tuple<typename Lubby2<DisplacementDim>::KelvinVector,
                           std::unique_ptr<typename MechanicsBase<
                               DisplacementDim>::MaterialStateVariables>,
                           typename Lubby2<DisplacementDim>::KelvinMatrix>>
Lubby2<DisplacementDim>::integrateStress(
    double const t,
    ProcessLib::SpatialPosition const& x,
    double const dt,
    KelvinVector const& /*eps_prev*/,
    KelvinVector const& eps,
    KelvinVector const& /*sigma_prev*/,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
        material_state_variables)
{
    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;

    assert(dynamic_cast<MaterialStateVariables const*>(
               &material_state_variables) != nullptr);
    MaterialStateVariables state(
        static_cast<MaterialStateVariables const&>(material_state_variables));
    state.setInitialConditions();

    auto local_lubby2_properties =
        detail::LocalLubby2Properties<DisplacementDim>{t, x, _mp};

    // calculation of deviatoric parts
    auto const& P_dev = Invariants::deviatoric_projection;
    KelvinVector const epsd_i = P_dev * eps;

    // initial guess as elastic predictor
    KelvinVector sigd_j = 2.0 * (epsd_i - state.eps_M_t - state.eps_K_t);

    // Calculate effective stress and update material properties
    double sig_eff = Invariants::equivalentStress(sigd_j);
    local_lubby2_properties.update(sig_eff);

    using LocalJacobianMatrix =
        Eigen::Matrix<double, KelvinVectorSize * 3, KelvinVectorSize * 3,
                      Eigen::RowMajor>;

    // Linear solver for the newton loop is required after the loop with the
    // same matrix. This saves one decomposition.
    Eigen::FullPivLU<LocalJacobianMatrix> linear_solver;

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

    LocalJacobianMatrix K_loc;
    {  // Local Newton solver
        using LocalResidualVector =
            Eigen::Matrix<double, KelvinVectorSize * 3, 1>;

        auto const update_residual = [&](LocalResidualVector& residual) {
            calculateResidualBurgers(
                dt, epsd_i, sigd_j, state.eps_K_j, state.eps_K_t, state.eps_M_j,
                state.eps_M_t, residual, local_lubby2_properties);
        };

        auto const update_jacobian = [&](LocalJacobianMatrix& jacobian) {
            calculateJacobianBurgers(
                t, x, dt, jacobian, sig_eff, sigd_j, state.eps_K_j,
                local_lubby2_properties);  // for solution dependent Jacobians
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
            local_lubby2_properties.update(sig_eff);
        };

        auto newton_solver = NumLib::NewtonRaphson<
            decltype(linear_solver), LocalJacobianMatrix,
            decltype(update_jacobian), LocalResidualVector,
            decltype(update_residual), decltype(update_solution)>(
            linear_solver, update_jacobian, update_residual, update_solution,
            _nonlinear_solver_parameters);

        auto const success_iterations = newton_solver.solve(K_loc);

        if (!success_iterations)
            return {};

        // If the Newton loop didn't run, the linear solver will not be
        // initialized.
        // This happens usually for the first iteration of the first timestep.
        if (*success_iterations == 0)
            linear_solver.compute(K_loc);
    }

    KelvinMatrix C =
        tangentStiffnessA<DisplacementDim>(local_lubby2_properties.GM0,
                                           local_lubby2_properties.KM0,
                                           linear_solver);

    // Hydrostatic part for the stress and the tangent.
    double const eps_i_trace = Invariants::trace(eps);
    return {
        {local_lubby2_properties.GM0 * sigd_j +
             local_lubby2_properties.KM0 * eps_i_trace * Invariants::identity2,
         std::unique_ptr<
             typename MechanicsBase<DisplacementDim>::MaterialStateVariables>{
             new MaterialStateVariables{state}},
         C}};
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
    detail::LocalLubby2Properties<DisplacementDim> const& properties)
{
    // calculate stress residual
    res.template block<KelvinVectorSize, 1>(0, 0).noalias() =
        stress_curr - 2. * (strain_curr - strain_Kel_curr - strain_Max_curr);

    // calculate Kelvin strain residual
    res.template block<KelvinVectorSize, 1>(KelvinVectorSize, 0).noalias() =
        1. / dt * (strain_Kel_curr - strain_Kel_t) -
        1. / (2. * properties.etaK) * (properties.GM0 * stress_curr -
                                       2. * properties.GK * strain_Kel_curr);

    // calculate Maxwell strain residual
    res.template block<KelvinVectorSize, 1>(2 * KelvinVectorSize, 0).noalias() =
        1. / dt * (strain_Max_curr - strain_Max_t) -
        0.5 * properties.GM0 / properties.etaM * stress_curr;
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
    detail::LocalLubby2Properties<DisplacementDim> const& properties)
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
        .noalias() =
        -0.5 * properties.GM0 / properties.etaK * KelvinMatrix::Identity();
    if (s_eff > 0.)
    {
        KelvinVector const eps_K_aid =
            1. / (properties.etaK * properties.etaK) *
            (properties.GM0 * sig_i - 2. * properties.GK * eps_K_i);

        KelvinVector const dG_K = 1.5 * _mp.mK(t, x)[0] * properties.GK *
                                  properties.GM0 / s_eff * sig_i;
        KelvinVector const dmu_vK = 1.5 * _mp.mvK(t, x)[0] * properties.GM0 *
                                    properties.etaK / s_eff * sig_i;
        Jac.template block<KelvinVectorSize, KelvinVectorSize>(KelvinVectorSize,
                                                               0)
            .noalias() += 0.5 * eps_K_aid * dmu_vK.transpose() +
                          1. / properties.etaK * eps_K_i * dG_K.transpose();
    }

    // build G_22
    Jac.template block<KelvinVectorSize, KelvinVectorSize>(KelvinVectorSize,
                                                           KelvinVectorSize)
        .diagonal()
        .setConstant(1. / dt + properties.GK / properties.etaK);

    // nothing to do for G_23

    // build G_31
    Jac.template block<KelvinVectorSize, KelvinVectorSize>(2 * KelvinVectorSize,
                                                           0)
        .noalias() =
        -0.5 * properties.GM0 / properties.etaM * KelvinMatrix::Identity();
    if (s_eff > 0.)
    {
        KelvinVector const dmu_vM = 1.5 * _mp.mvM(t, x)[0] * properties.GM0 *
                                    properties.etaM / s_eff * sig_i;
        Jac.template block<KelvinVectorSize, KelvinVectorSize>(
               2 * KelvinVectorSize, 0)
            .noalias() += 0.5 * properties.GM0 /
                          (properties.etaM * properties.etaM) * sig_i *
                          dmu_vM.transpose();
    }

    // nothing to do for G_32

    // build G_33
    Jac.template block<KelvinVectorSize, KelvinVectorSize>(2 * KelvinVectorSize,
                                                           2 * KelvinVectorSize)
        .diagonal()
        .setConstant(1. / dt);
}

}  // namespace Lubby2
}  // namespace Solids
}  // namespace MaterialLib
