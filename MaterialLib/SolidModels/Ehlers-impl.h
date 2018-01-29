/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**
 * Common convenitions for naming:
 * x_D              - deviatoric part of tensor x
 * x_V              - volumetric part of tensor x
 * x_p              - a variable related to plastic potential
 * x_prev           - value of x in previous time step
 *
 * Variables used in the code:
 * eps_D            - deviatoric strain
 * eps_p_D_dot      - deviatoric increment of plastic strain
 * eps_p_eff_dot    - increment of effective plastic strain
 * eps_p_V_dot      - volumetric increment of plastic strain
 * sigma_D_inverse_D - deviatoric part of sigma_D_inverse
 *
 * derivation of the flow rule
 * theta            - J3 / J2^(3 / 2) from yield function
 * dtheta_dsigma    - derivative of theta
 * sqrtPhi          - square root of Phi from plastic potential
 * flow_D           - deviatoric part of flow
 * flow_V           - volumetric part of flow
 * lambda_flow_D    - deviatoric increment of plastic strain
 *
 */
#pragma once

#include <boost/math/special_functions/pow.hpp>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace MaterialLib
{
namespace Solids
{
namespace Ehlers
{
/// Special product of \c v with itself: \f$v \odot v\f$.
/// The tensor \c v is given in Kelvin mapping.
/// \note Implementation only for 2 and 3 dimensions.
/// \attention Pay attention to the sign of the result, which normally would be
/// negative, but the returned value is not negated. This has to do with \f$
/// d(A^{-1})/dA = -A^{-1} \odot A^{-1} \f$.
template <int DisplacementDim>
MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> sOdotS(
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& v);

template <int DisplacementDim>
struct PhysicalStressWithInvariants final
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

    explicit PhysicalStressWithInvariants(KelvinVector const& stress)
        : value{stress},
          D{Invariants::deviatoric_projection * stress},
          I_1{Invariants::trace(stress)},
          J_2{Invariants::J2(D)},
          J_3{Invariants::J3(D)}
    {
    }

    KelvinVector value;
    KelvinVector D;
    double I_1;
    double J_2;
    double J_3;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

/// Holds powers of 1 + gamma_p*theta to base 0, m_p, and m_p-1.
struct OnePlusGamma_pTheta final
{
    OnePlusGamma_pTheta(double const gamma_p, double const theta,
                        double const m_p)
        : value{1 + gamma_p * theta},
          pow_m_p{std::pow(value, m_p)},
          pow_m_p1{pow_m_p / value}
    {
    }

    double const value;
    double const pow_m_p;
    double const pow_m_p1;
};

template <int DisplacementDim>
double plasticFlowVolumetricPart(
    PhysicalStressWithInvariants<DisplacementDim> const& s,
    double const sqrtPhi, double const alpha_p, double const beta_p,
    double const delta_p, double const epsilon_p)
{
    return 3 * (alpha_p * s.I_1 +
                4 * boost::math::pow<2>(delta_p) * boost::math::pow<3>(s.I_1)) /
               (2 * sqrtPhi) +
           3 * beta_p + 6 * epsilon_p * s.I_1;
}

template <int DisplacementDim>
typename SolidEhlers<DisplacementDim>::KelvinVector plasticFlowDeviatoricPart(
    PhysicalStressWithInvariants<DisplacementDim> const& s,
    OnePlusGamma_pTheta const& one_gt, double const sqrtPhi,
    typename SolidEhlers<DisplacementDim>::KelvinVector const& dtheta_dsigma,
    double const gamma_p, double const m_p)
{
    return (one_gt.pow_m_p *
            (s.D + s.J_2 * m_p * gamma_p * dtheta_dsigma / one_gt.value)) /
           (2 * sqrtPhi);
}

template <int DisplacementDim>
double yieldFunction(MaterialProperties const& mp,
                     PhysicalStressWithInvariants<DisplacementDim> const& s,
                     double const k)
{
    double const I_1_squared = boost::math::pow<2>(s.I_1);
    assert(s.J_2 != 0);

    return std::sqrt(
               s.J_2 *
                   std::pow(1 + mp.gamma * s.J_3 / (s.J_2 * std::sqrt(s.J_2)),
                            mp.m) +
               mp.alpha / 2. * I_1_squared +
               boost::math::pow<2>(mp.delta) *
                   boost::math::pow<2>(I_1_squared)) +
           mp.beta * s.I_1 + mp.epsilon * I_1_squared - k;
}

template <int DisplacementDim>
typename SolidEhlers<DisplacementDim>::ResidualVectorType
calculatePlasticResidual(
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps_D,
    double const eps_V,
    PhysicalStressWithInvariants<DisplacementDim> const& s,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps_p_D,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps_p_D_dot,
    double const eps_p_V,
    double const eps_p_V_dot,
    double const eps_p_eff_dot,
    double const lambda,
    double const k,
    MaterialProperties const& mp)
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

    auto const& P_dev = Invariants::deviatoric_projection;
    auto const& identity2 = Invariants::identity2;

    double const theta = s.J_3 / (s.J_2 * std::sqrt(s.J_2));

    typename SolidEhlers<DisplacementDim>::ResidualVectorType residual;
    // calculate stress residual
    residual.template segment<KelvinVectorSize>(0).noalias() =
        s.value / mp.G - 2 * (eps_D - eps_p_D) -
        mp.K / mp.G * (eps_V - eps_p_V) * identity2;

    // deviatoric plastic strain
    KelvinVector const sigma_D_inverse_D =
        P_dev * MathLib::KelvinVector::inverse(s.D);
    KelvinVector const dtheta_dsigma =
        theta * sigma_D_inverse_D - 3. / 2. * theta / s.J_2 * s.D;

    OnePlusGamma_pTheta const one_gt{mp.gamma_p, theta, mp.m_p};
    double const sqrtPhi = std::sqrt(
        s.J_2 * one_gt.pow_m_p + mp.alpha_p / 2. * boost::math::pow<2>(s.I_1) +
        boost::math::pow<2>(mp.delta_p) * boost::math::pow<4>(s.I_1));
    KelvinVector const flow_D = plasticFlowDeviatoricPart(
        s, one_gt, sqrtPhi, dtheta_dsigma, mp.gamma_p, mp.m_p);
    KelvinVector const lambda_flow_D = lambda * flow_D;

    residual.template segment<KelvinVectorSize>(KelvinVectorSize).noalias() =
        eps_p_D_dot - lambda_flow_D;

    // plastic volume strain
    {
        double const flow_V = plasticFlowVolumetricPart<DisplacementDim>(
            s, sqrtPhi, mp.alpha_p, mp.beta_p, mp.delta_p, mp.epsilon_p);
        residual(2 * KelvinVectorSize, 0) = eps_p_V_dot - lambda * flow_V;
    }

    // evolution of plastic equivalent strain
    residual(2 * KelvinVectorSize + 1) =
        eps_p_eff_dot -
        std::sqrt(2. / 3. * lambda_flow_D.transpose() * lambda_flow_D);

    // yield function (for plastic multiplier)
    residual(2 * KelvinVectorSize + 2) = yieldFunction(mp, s, k) / mp.G;
    return residual;
}

template <int DisplacementDim>
typename SolidEhlers<DisplacementDim>::JacobianMatrix calculatePlasticJacobian(
    double const dt,
    PhysicalStressWithInvariants<DisplacementDim> const& s,
    double const lambda,
    MaterialProperties const& mp)
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    auto const& P_dev = Invariants::deviatoric_projection;
    auto const& identity2 = Invariants::identity2;

    double const theta = s.J_3 / (s.J_2 * std::sqrt(s.J_2));
    OnePlusGamma_pTheta const one_gt{mp.gamma_p, theta, mp.m_p};

    // inverse of deviatoric stress tensor
    if (Invariants::determinant(s.D) == 0)
    {
        OGS_FATAL("Determinant is zero. Matrix is non-invertable.");
    }
    // inverse of sigma_D
    KelvinVector const sigma_D_inverse = MathLib::KelvinVector::inverse(s.D);
    KelvinVector const sigma_D_inverse_D = P_dev * sigma_D_inverse;

    KelvinVector const dtheta_dsigma =
        theta * sigma_D_inverse_D - 3. / 2. * theta / s.J_2 * s.D;

    // deviatoric flow
    double const sqrtPhi = std::sqrt(
        s.J_2 * one_gt.pow_m_p + mp.alpha_p / 2. * boost::math::pow<2>(s.I_1) +
        boost::math::pow<2>(mp.delta_p) * boost::math::pow<4>(s.I_1));
    KelvinVector const flow_D = plasticFlowDeviatoricPart(
        s, one_gt, sqrtPhi, dtheta_dsigma, mp.gamma_p, mp.m_p);
    KelvinVector const lambda_flow_D = lambda * flow_D;

    typename SolidEhlers<DisplacementDim>::JacobianMatrix jacobian =
        SolidEhlers<DisplacementDim>::JacobianMatrix::Zero();

    // G_11
    jacobian.template block<KelvinVectorSize, KelvinVectorSize>(0, 0)
        .noalias() = KelvinMatrix::Identity();

    // G_12
    jacobian
        .template block<KelvinVectorSize, KelvinVectorSize>(0, KelvinVectorSize)
        .noalias() = 2 * KelvinMatrix::Identity();

    // G_13
    jacobian.template block<KelvinVectorSize, 1>(0, 2 * KelvinVectorSize)
        .noalias() = mp.K / mp.G * identity2;

    // G_14 and G_15 are zero

    // G_21 -- derivative of deviatoric flow

    double const gm_p = mp.gamma_p * mp.m_p;
    // intermediate variable for derivative of deviatoric flow
    KelvinVector const M0 = s.J_2 / one_gt.value * dtheta_dsigma;
    // derivative of Phi w.r.t. sigma
    KelvinVector const dPhi_dsigma =
        one_gt.pow_m_p * (s.D + gm_p * M0) +
        (mp.alpha_p * s.I_1 +
         4 * boost::math::pow<2>(mp.delta_p) * boost::math::pow<3>(s.I_1)) *
            identity2;

    // intermediate variable for derivative of deviatoric flow
    KelvinMatrix const M1 =
        one_gt.pow_m_p *
        (s.D * dPhi_dsigma.transpose() + gm_p * M0 * dPhi_dsigma.transpose());
    // intermediate variable for derivative of deviatoric flow
    KelvinMatrix const M2 =
        one_gt.pow_m_p * (P_dev + s.D * gm_p * M0.transpose());
    // second derivative of theta
    KelvinMatrix const d2theta_dsigma2 =
        theta * P_dev * sOdotS<DisplacementDim>(sigma_D_inverse) * P_dev +
        sigma_D_inverse_D * dtheta_dsigma.transpose() -
        3. / 2. * theta / s.J_2 * P_dev -
        3. / 2. * dtheta_dsigma / s.J_2 * s.D.transpose() +
        3. / 2. * theta / boost::math::pow<2>(s.J_2) * s.D * s.D.transpose();

    // intermediate variable for derivative of deviatoric flow
    KelvinMatrix const M3 =
        gm_p * one_gt.pow_m_p1 *
        ((s.D + (gm_p - mp.gamma_p) * M0) * dtheta_dsigma.transpose() +
         s.J_2 * d2theta_dsigma2);

    // derivative of flow_D w.r.t. sigma
    KelvinMatrix const dflow_D_dsigma =
        (-M1 / (4 * boost::math::pow<3>(sqrtPhi)) + (M2 + M3) / (2 * sqrtPhi)) *
        mp.G;
    jacobian
        .template block<KelvinVectorSize, KelvinVectorSize>(KelvinVectorSize, 0)
        .noalias() = -lambda * dflow_D_dsigma;

    // G_22
    jacobian
        .template block<KelvinVectorSize, KelvinVectorSize>(KelvinVectorSize,
                                                            KelvinVectorSize)
        .noalias() = KelvinMatrix::Identity() / dt;

    // G_23 and G_24 are zero

    // G_25
    jacobian
        .template block<KelvinVectorSize, 1>(KelvinVectorSize,
                                             2 * KelvinVectorSize + 2)
        .noalias() = -flow_D;

    // G_31
    {
        // derivative of flow_V w.r.t. sigma
        KelvinVector const dflow_V_dsigma =
            3 * mp.G *
            (-(mp.alpha_p * s.I_1 +
               4 * boost::math::pow<2>(mp.delta_p) *
                   boost::math::pow<3>(s.I_1)) /
                 (4 * boost::math::pow<3>(sqrtPhi)) * dPhi_dsigma +
             (mp.alpha_p * identity2 +
              12 * boost::math::pow<2>(mp.delta_p * s.I_1) * identity2) /
                 (2 * sqrtPhi) +
             2 * mp.epsilon_p * identity2);

        jacobian.template block<1, KelvinVectorSize>(2 * KelvinVectorSize, 0)
            .noalias() = -lambda * dflow_V_dsigma.transpose();
    }

    // G_32 is zero

    // G_33
    jacobian(2 * KelvinVectorSize, 2 * KelvinVectorSize) = 1. / dt;

    // G_34 is zero

    // G_35
    {
        double const flow_V = plasticFlowVolumetricPart<DisplacementDim>(
            s, sqrtPhi, mp.alpha_p, mp.beta_p, mp.delta_p, mp.epsilon_p);
        jacobian(2 * KelvinVectorSize, 2 * KelvinVectorSize + 2) = -flow_V;
    }

    // increment of effectiv plastic strain
    double const eff_flow =
        std::sqrt(2. / 3. * lambda_flow_D.transpose() * lambda_flow_D);

    if (eff_flow > 0)
    {
        // intermediate variable for derivative of plastic jacobian
        KelvinVector const eff_flow23_lambda_flow_D =
            -2 / 3. / eff_flow * lambda_flow_D;
        // G_41
        jacobian
            .template block<1, KelvinVectorSize>(2 * KelvinVectorSize + 1, 0)
            .noalias() = lambda * dflow_D_dsigma * eff_flow23_lambda_flow_D;
        // G_45
        jacobian(2 * KelvinVectorSize + 1, 2 * KelvinVectorSize + 2) =
            eff_flow23_lambda_flow_D.transpose() * flow_D;
    }

    // G_42 and G_43 are zero

    // G_44
    jacobian(2 * KelvinVectorSize + 1, 2 * KelvinVectorSize + 1) = 1. / dt;

    // G_51
    {
        double const one_gt_pow_m = std::pow(one_gt.value, mp.m);
        double const gm = mp.gamma * mp.m;
        // derivative of yield function w.r.t. sigma
        KelvinVector const dF_dsigma =
            mp.G * (one_gt_pow_m * (s.D + gm * M0) +
                    (mp.alpha * s.I_1 +
                     4 * boost::math::pow<2>(mp.delta) *
                         boost::math::pow<3>(s.I_1)) *
                        identity2) /
                (2. * sqrtPhi) +
            mp.G * (mp.beta + 2 * mp.epsilon_p * s.I_1) * identity2;

        jacobian
            .template block<1, KelvinVectorSize>(2 * KelvinVectorSize + 2, 0)
            .noalias() = dF_dsigma.transpose() / mp.G;
    }

    // G_54
    jacobian(2 * KelvinVectorSize + 2, 2 * KelvinVectorSize + 1) =
        -mp.kappa * mp.hardening_coefficient / mp.G;

    // G_52, G_53, G_55 are zero
    return jacobian;
}

/// Calculates the derivative of the residuals with respect to total
/// strain. Implementation fully implicit only.
template <int DisplacementDim>
MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> calculateDResidualDEps(
    double const K, double const G)
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    auto const& P_dev = Invariants::deviatoric_projection;
    auto const& P_sph = Invariants::spherical_projection;
    auto const& I =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>::Identity();

    return -2. * I * P_dev - 3. * K / G * I * P_sph;
}

inline double calculateIsotropicHardening(double const kappa,
                                          double const hardening_coefficient,
                                          double const eps_p_eff)
{
    return kappa * (1. + eps_p_eff * hardening_coefficient);
}

template <int DisplacementDim>
typename SolidEhlers<DisplacementDim>::KelvinVector predict_sigma(
    double const G, double const K,
    typename SolidEhlers<DisplacementDim>::KelvinVector const& sigma_prev,
    typename SolidEhlers<DisplacementDim>::KelvinVector const& eps,
    typename SolidEhlers<DisplacementDim>::KelvinVector const& eps_prev,
    double const eps_V)
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    auto const& P_dev = Invariants::deviatoric_projection;

    // dimensionless initial hydrostatic pressure
    double const pressure_prev = Invariants::trace(sigma_prev) / (-3. * G);
    // initial strain invariant
    double const e_prev = Invariants::trace(eps_prev);
    // dimensioness hydrostatic stress increment
    double const pressure = pressure_prev - K / G * (eps_V - e_prev);
    // dimensionless deviatoric initial stress
    typename SolidEhlers<DisplacementDim>::KelvinVector const sigma_D_prev =
        P_dev * sigma_prev / G;
    // dimensionless deviatoric stress
    typename SolidEhlers<DisplacementDim>::KelvinVector const sigma_D =
        sigma_D_prev + 2 * P_dev * (eps - eps_prev);
    return sigma_D - pressure * Invariants::identity2;
}

/// Split the agglomerated solution vector in separate items. The arrangement
/// must be the same as in the newton() function.
template <typename ResidualVector, typename KelvinVector>
std::tuple<KelvinVector, PlasticStrain<KelvinVector>, double>
splitSolutionVector(ResidualVector const& solution)
{
    static auto const size = KelvinVector::SizeAtCompileTime;
    return std::forward_as_tuple(
        solution.template segment<size>(size * 0),
        PlasticStrain<KelvinVector>{solution.template segment<size>(size * 1),
                                    solution[size * 2], solution[size * 2 + 1]},
        solution[size * 2 + 2]);
}

template <int DisplacementDim>
boost::optional<std::tuple<typename SolidEhlers<DisplacementDim>::KelvinVector,
                           std::unique_ptr<typename MechanicsBase<
                               DisplacementDim>::MaterialStateVariables>,
                           typename SolidEhlers<DisplacementDim>::KelvinMatrix>>
SolidEhlers<DisplacementDim>::integrateStress(
    double const t,
    ProcessLib::SpatialPosition const& x,
    double const dt,
    KelvinVector const& eps_prev,
    KelvinVector const& eps,
    KelvinVector const& sigma_prev,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
        material_state_variables)
{
    assert(dynamic_cast<StateVariables<DisplacementDim> const*>(
               &material_state_variables) != nullptr);

    StateVariables<DisplacementDim> state =
        static_cast<StateVariables<DisplacementDim> const&>(
            material_state_variables);
    state.setInitialConditions();

    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    // volumetric strain
    double const eps_V = Invariants::trace(eps);

    auto const& P_dev = Invariants::deviatoric_projection;
    // deviatoric strain
    KelvinVector const eps_D = P_dev * eps;

    // do the evaluation once per function call.
    MaterialProperties const mp(t, x, _mp);

    KelvinVector sigma = predict_sigma<DisplacementDim>(mp.G, mp.K, sigma_prev,
                                                        eps, eps_prev, eps_V);

    KelvinMatrix tangentStiffness;

    PhysicalStressWithInvariants<DisplacementDim> s{mp.G * sigma};
    // Quit early if sigma is zero (nothing to do) or if we are still in elastic
    // zone.
    if ((sigma.squaredNorm() == 0 ||
         yieldFunction(
             mp, s,
             calculateIsotropicHardening(mp.kappa, mp.hardening_coefficient,
                                         state.eps_p.eff)) < 0))
    {
        tangentStiffness.setZero();
        tangentStiffness.template topLeftCorner<3, 3>().setConstant(
            mp.K - 2. / 3 * mp.G);
        tangentStiffness.noalias() += 2 * mp.G * KelvinMatrix::Identity();
    }
    else
    {
        // Linear solver for the newton loop is required after the loop with the
        // same matrix. This saves one decomposition.
        Eigen::FullPivLU<Eigen::Matrix<double, JacobianResidualSize,
                                       JacobianResidualSize, Eigen::RowMajor>>
            linear_solver;

        {
            static int const KelvinVectorSize =
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value;
            using KelvinVector =
                MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
            using ResidualVectorType =
                Eigen::Matrix<double, JacobianResidualSize, 1>;
            using JacobianMatrix =
                Eigen::Matrix<double, JacobianResidualSize,
                              JacobianResidualSize, Eigen::RowMajor>;

            JacobianMatrix jacobian;

            // Agglomerated solution vector construction.  It is later split
            // into individual parts by splitSolutionVector().
            ResidualVectorType solution;
            solution << sigma, state.eps_p.D, state.eps_p.V, state.eps_p.eff, 0;

            auto const update_residual = [&](ResidualVectorType& residual) {

                auto const& eps_p_D =
                    solution.template segment<KelvinVectorSize>(
                        KelvinVectorSize);
                KelvinVector const eps_p_D_dot =
                    (eps_p_D - state.eps_p_prev.D) / dt;

                double const& eps_p_V = solution[KelvinVectorSize * 2];
                double const eps_p_V_dot =
                    (eps_p_V - state.eps_p_prev.V) / dt;

                double const& eps_p_eff = solution[KelvinVectorSize * 2 + 1];
                double const eps_p_eff_dot =
                    (eps_p_eff - state.eps_p_prev.eff) / dt;

                double const k_hardening = calculateIsotropicHardening(
                    mp.kappa, mp.hardening_coefficient,
                    solution[KelvinVectorSize * 2 + 1]);
                residual = calculatePlasticResidual<DisplacementDim>(
                    eps_D, eps_V, s,
                    solution.template segment<KelvinVectorSize>(
                        KelvinVectorSize),
                    eps_p_D_dot, solution[KelvinVectorSize * 2], eps_p_V_dot,
                    eps_p_eff_dot, solution[KelvinVectorSize * 2 + 2],
                    k_hardening, mp);
            };

            auto const update_jacobian = [&](JacobianMatrix& jacobian) {
                jacobian = calculatePlasticJacobian<DisplacementDim>(
                    dt, s, solution[KelvinVectorSize * 2 + 2], mp);
            };

            auto const update_solution =
                [&](ResidualVectorType const& increment) {
                    solution += increment;
                    s = PhysicalStressWithInvariants<DisplacementDim>{
                        mp.G * solution.template segment<KelvinVectorSize>(0)};
                };

            auto newton_solver = NumLib::NewtonRaphson<
                decltype(linear_solver), JacobianMatrix,
                decltype(update_jacobian), ResidualVectorType,
                decltype(update_residual), decltype(update_solution)>(
                linear_solver, update_jacobian, update_residual,
                update_solution, _nonlinear_solver_parameters);

            auto const success_iterations = newton_solver.solve(jacobian);

            if (!success_iterations)
                return {};

            // If the Newton loop didn't run, the linear solver will not be
            // initialized.
            // This happens usually for the first iteration of the first
            // timestep.
            if (*success_iterations == 0)
                linear_solver.compute(jacobian);

            std::tie(sigma, state.eps_p, std::ignore) =
                splitSolutionVector<ResidualVectorType, KelvinVector>(solution);
        }

        // Calculate residual derivative w.r.t. strain
        Eigen::Matrix<double, JacobianResidualSize, KelvinVectorSize,
                      Eigen::RowMajor>
            dresidual_deps =
                Eigen::Matrix<double, JacobianResidualSize, KelvinVectorSize,
                              Eigen::RowMajor>::Zero();
        dresidual_deps.template block<KelvinVectorSize, KelvinVectorSize>(0, 0)
            .noalias() = calculateDResidualDEps<DisplacementDim>(mp.K, mp.G);

        tangentStiffness =
            mp.G *
            linear_solver.solve(-dresidual_deps)
                .template block<KelvinVectorSize, KelvinVectorSize>(0, 0);
    }

    KelvinVector sigma_final = mp.G * sigma;

    return {std::make_tuple(
        sigma_final,
        std::unique_ptr<
            typename MechanicsBase<DisplacementDim>::MaterialStateVariables>{
            new StateVariables<DisplacementDim>{
                static_cast<StateVariables<DisplacementDim> const&>(state)}},
        tangentStiffness)};
}

template <int DisplacementDim>
std::vector<typename MechanicsBase<DisplacementDim>::InternalVariable>
SolidEhlers<DisplacementDim>::getInternalVariables() const
{
    return {
        {"damage.kappa_d", 1,
         [](typename MechanicsBase<
                DisplacementDim>::MaterialStateVariables const& state,
            std::vector<double>& cache) -> std::vector<double> const& {
             assert(dynamic_cast<StateVariables<DisplacementDim> const*>(
                        &state) != nullptr);
             auto const& ehlers_state =
                 static_cast<StateVariables<DisplacementDim> const&>(state);

             cache.resize(1);
             cache.front() = ehlers_state.damage.kappa_d();
             return cache;
         }},
        {"damage.value", 1,
         [](typename MechanicsBase<
                DisplacementDim>::MaterialStateVariables const& state,
            std::vector<double>& cache) -> std::vector<double> const& {
             assert(dynamic_cast<StateVariables<DisplacementDim> const*>(
                        &state) != nullptr);
             auto const& ehlers_state =
                 static_cast<StateVariables<DisplacementDim> const&>(state);

             cache.resize(1);
             cache.front() = ehlers_state.damage.value();
             return cache;
         }},
        {"eps_p.D", KelvinVector::RowsAtCompileTime,
         [](typename MechanicsBase<
                DisplacementDim>::MaterialStateVariables const& state,
            std::vector<double>& cache) -> std::vector<double> const& {
             assert(dynamic_cast<StateVariables<DisplacementDim> const*>(
                        &state) != nullptr);
             auto const& ehlers_state =
                 static_cast<StateVariables<DisplacementDim> const&>(state);

             cache.resize(KelvinVector::RowsAtCompileTime);
             MathLib::toVector<KelvinVector>(cache,
                                             KelvinVector::RowsAtCompileTime) =
                 MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                     ehlers_state.eps_p.D);

             return cache;
         }},
        {"eps_p.V", 1,
         [](typename MechanicsBase<
                DisplacementDim>::MaterialStateVariables const& state,
            std::vector<double>& cache) -> std::vector<double> const& {
             assert(dynamic_cast<StateVariables<DisplacementDim> const*>(
                        &state) != nullptr);
             auto const& ehlers_state =
                 static_cast<StateVariables<DisplacementDim> const&>(state);

             cache.resize(1);
             cache.front() = ehlers_state.eps_p.V;
             return cache;
         }},
        {"eps_p.eff", 1,
         [](typename MechanicsBase<
                DisplacementDim>::MaterialStateVariables const& state,
            std::vector<double>& cache) -> std::vector<double> const& {
             assert(dynamic_cast<StateVariables<DisplacementDim> const*>(
                        &state) != nullptr);
             auto const& ehlers_state =
                 static_cast<StateVariables<DisplacementDim> const&>(state);

             cache.resize(1);
             cache.front() = ehlers_state.eps_p.eff;
             return cache;
         }}};
}

}  // namespace Ehlers
}  // namespace Solids
}  // namespace MaterialLib
