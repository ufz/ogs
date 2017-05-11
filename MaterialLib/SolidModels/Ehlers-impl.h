/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
#include <logog/include/logog.hpp>
#include "MaterialLib/SolidModels/KelvinVector.h"

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
ProcessLib::KelvinMatrixType<DisplacementDim> sOdotS(
    ProcessLib::KelvinVectorType<DisplacementDim> const& v);

template <int DisplacementDim>
struct PhysicalStressWithInvariants final
{
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;

    explicit PhysicalStressWithInvariants(KelvinVector const& stress)
        : value{stress},
          D{Invariants::deviatoric_projection * stress},
          I_1{Invariants::trace(stress)},
          J_2{Invariants::J2(D)},
          J_3{Invariants::J3(D)}
    {
    }

    PhysicalStressWithInvariants(PhysicalStressWithInvariants const&) = default;
    PhysicalStressWithInvariants& operator=(
        PhysicalStressWithInvariants const&) = default;
#if defined(_MSC_VER) && (_MSC_VER >= 1900)
    PhysicalStressWithInvariants(PhysicalStressWithInvariants&&) = default;
    PhysicalStressWithInvariants& operator=(PhysicalStressWithInvariants&&) =
        default;
#endif  // _MSC_VER

    KelvinVector value;
    KelvinVector D;
    double I_1;
    double J_2;
    double J_3;
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
double yieldFunction(PhysicalStressWithInvariants<DisplacementDim> const& s,
                     MaterialProperties const& mp, double const t,
                     ProcessLib::SpatialPosition const& x, double const k)
{
    double const alpha = mp.alpha(t, x)[0];
    double const beta = mp.beta(t, x)[0];
    double const delta = mp.delta(t, x)[0];
    double const epsilon = mp.epsilon(t, x)[0];
    double const gamma = mp.gamma(t, x)[0];
    double const m = mp.m(t, x)[0];

    double const I_1_squared = boost::math::pow<2>(s.I_1);
    assert(s.J_2 != 0);

    return std::sqrt(
               s.J_2 * std::pow(1 +
                                    gamma * s.J_3 /
                                        boost::math::pow<3>(std::sqrt(s.J_2)),
                                m) +
               alpha / 2. * I_1_squared +
               boost::math::pow<2>(delta) * boost::math::pow<2>(I_1_squared)) +
           beta * s.I_1 + epsilon * I_1_squared - k;
}

template <int DisplacementDim>
void calculatePlasticResidual(
    double const t,
    ProcessLib::SpatialPosition const& x,
    ProcessLib::KelvinVectorType<DisplacementDim> const& eps_D,
    double const eps_V,
    PhysicalStressWithInvariants<DisplacementDim> const& s,
    ProcessLib::KelvinVectorType<DisplacementDim> const& eps_p_D,
    ProcessLib::KelvinVectorType<DisplacementDim> const& eps_p_D_dot,
    double const eps_p_V,
    double const eps_p_V_dot,
    double const eps_p_eff_dot,
    double const lambda,
    double const k,
    MaterialProperties const& _mp,
    typename SolidEhlers<DisplacementDim>::ResidualVectorType& residual)
{
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;

    auto const& P_dev = Invariants::deviatoric_projection;
    auto const& identity2 = Invariants::identity2;

    double const G = _mp.G(t, x)[0];
    double const K = _mp.K(t, x)[0];

    double const theta = s.J_3 / boost::math::pow<3>(std::sqrt(s.J_2));

    // calculate stress residual
    residual.template segment<KelvinVectorSize>(0).noalias() =
        s.value / G - 2 * (eps_D - eps_p_D) -
        K / G * (eps_V - eps_p_V) * identity2;

    // deviatoric plastic strain
    double const alpha_p = _mp.alpha_p(t, x)[0];
    double const delta_p = _mp.delta_p(t, x)[0];
    double const gamma_p = _mp.gamma_p(t, x)[0];
    double const m_p = _mp.m_p(t, x)[0];
    KelvinVector const sigma_D_inverse_D =
        P_dev * MaterialLib::SolidModels::inverse(s.D);
    KelvinVector const dtheta_dsigma =
        theta * sigma_D_inverse_D - 3. / 2. * theta / s.J_2 * s.D;

    OnePlusGamma_pTheta const one_gt{gamma_p, theta, m_p};
    double const sqrtPhi = std::sqrt(
        s.J_2 * one_gt.pow_m_p + alpha_p / 2. * boost::math::pow<2>(s.I_1) +
        boost::math::pow<2>(delta_p) * boost::math::pow<4>(s.I_1));
    KelvinVector const flow_D = plasticFlowDeviatoricPart(
        s, one_gt, sqrtPhi, dtheta_dsigma, gamma_p, m_p);
    KelvinVector const lambda_flow_D = lambda * flow_D;

    residual.template segment<KelvinVectorSize>(KelvinVectorSize).noalias() =
        eps_p_D_dot - lambda_flow_D;

    // plastic volume strain
    {
        double const beta_p = _mp.beta_p(t, x)[0];
        double const epsilon_p = _mp.epsilon_p(t, x)[0];

        double const flow_V = plasticFlowVolumetricPart<DisplacementDim>(
            s, sqrtPhi, alpha_p, beta_p, delta_p, epsilon_p);
        residual(2 * KelvinVectorSize, 0) = eps_p_V_dot - lambda * flow_V;
    }

    // evolution of plastic equivalent strain
    residual(2 * KelvinVectorSize + 1) =
        eps_p_eff_dot -
        std::sqrt(2. / 3. * lambda_flow_D.transpose() * lambda_flow_D);

    // yield function (for plastic multiplier)
    residual(2 * KelvinVectorSize + 2) =
        yieldFunction<DisplacementDim>(s, _mp, t, x, k) / G;
}

template <int DisplacementDim>
void calculatePlasticJacobian(
    double const dt,
    double const t,
    ProcessLib::SpatialPosition const& x,
    typename SolidEhlers<DisplacementDim>::JacobianMatrix& jacobian,
    PhysicalStressWithInvariants<DisplacementDim> const& s,
    double const lambda,
    MaterialProperties const& _mp)
{
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

    auto const& P_dev = Invariants::deviatoric_projection;
    auto const& identity2 = Invariants::identity2;

    double const G = _mp.G(t, x)[0];
    double const K = _mp.K(t, x)[0];
    double const m = _mp.m(t, x)[0];
    double const alpha = _mp.alpha(t, x)[0];
    double const beta = _mp.beta(t, x)[0];
    double const gamma = _mp.gamma(t, x)[0];
    double const delta = _mp.delta(t, x)[0];

    double const alpha_p = _mp.alpha_p(t, x)[0];
    double const beta_p = _mp.beta_p(t, x)[0];
    double const gamma_p = _mp.gamma_p(t, x)[0];
    double const delta_p = _mp.delta_p(t, x)[0];
    double const epsilon_p = _mp.epsilon_p(t, x)[0];
    double const m_p = _mp.m_p(t, x)[0];

    double const theta = s.J_3 / boost::math::pow<3>(std::sqrt(s.J_2));
    OnePlusGamma_pTheta const one_gt{gamma_p, theta, m_p};

    // inverse of deviatoric stress tensor
    if (Invariants::determinant(s.D) == 0)
    {
        OGS_FATAL("Determinant is zero. Matrix is non-invertable.");
    }
    // inverse of sigma_D
    KelvinVector const sigma_D_inverse = MaterialLib::SolidModels::inverse(s.D);
    KelvinVector const sigma_D_inverse_D = P_dev * sigma_D_inverse;

    KelvinVector const dtheta_dsigma =
        theta * sigma_D_inverse_D - 3. / 2. * theta / s.J_2 * s.D;

    // deviatoric flow
    double const sqrtPhi = std::sqrt(
        s.J_2 * one_gt.pow_m_p + alpha_p / 2. * boost::math::pow<2>(s.I_1) +
        boost::math::pow<2>(delta_p) * boost::math::pow<4>(s.I_1));
    KelvinVector const flow_D = plasticFlowDeviatoricPart(
        s, one_gt, sqrtPhi, dtheta_dsigma, gamma_p, m_p);
    KelvinVector const lambda_flow_D = lambda * flow_D;

    jacobian.setZero();

    // G_11
    jacobian.template block<KelvinVectorSize, KelvinVectorSize>(0, 0)
        .noalias() = KelvinMatrix::Identity();

    // G_12
    jacobian
        .template block<KelvinVectorSize, KelvinVectorSize>(0, KelvinVectorSize)
        .noalias() = 2 * KelvinMatrix::Identity();

    // G_13
    jacobian.template block<KelvinVectorSize, 1>(0, 2 * KelvinVectorSize)
        .noalias() = K / G * identity2;

    // G_14 and G_15 are zero

    // G_21 -- derivative of deviatoric flow

    double const gm_p = gamma_p * m_p;
    // intermediate variable for derivative of deviatoric flow
    KelvinVector const M0 = s.J_2 / one_gt.value * dtheta_dsigma;
    // derivative of Phi w.r.t. sigma
    KelvinVector const dPhi_dsigma =
        one_gt.pow_m_p * (s.D + gm_p * M0) +
        (alpha_p * s.I_1 +
         4 * boost::math::pow<2>(delta_p) * boost::math::pow<3>(s.I_1)) *
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
        ((s.D + (gm_p - gamma_p) * M0) * dtheta_dsigma.transpose() +
         s.J_2 * d2theta_dsigma2);

    // derivative of flow_D w.r.t. sigma
    KelvinMatrix const dflow_D_dsigma =
        (-M1 / (4 * boost::math::pow<3>(sqrtPhi)) + (M2 + M3) / (2 * sqrtPhi)) *
        G;
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
            3 * G *
            (-(alpha_p * s.I_1 +
               4 * boost::math::pow<2>(delta_p) * boost::math::pow<3>(s.I_1)) /
                 (4 * boost::math::pow<3>(sqrtPhi)) * dPhi_dsigma +
             (alpha_p * identity2 +
              12 * boost::math::pow<2>(delta_p * s.I_1) * identity2) /
                 (2 * sqrtPhi) +
             2 * epsilon_p * identity2);

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
            s, sqrtPhi, alpha_p, beta_p, delta_p, epsilon_p);
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
        double const one_gt_pow_m = std::pow(one_gt.value, m);
        double const gm = gamma * m;
        // derivative of yield function w.r.t. sigma
        KelvinVector const dF_dsigma =
            G * (one_gt_pow_m * (s.D + gm * M0) +
                 (alpha * s.I_1 +
                  4 * boost::math::pow<2>(delta) * boost::math::pow<3>(s.I_1)) *
                     identity2) /
                (2. * sqrtPhi) +
            G * (beta + 2 * epsilon_p * s.I_1) * identity2;

        jacobian
            .template block<1, KelvinVectorSize>(2 * KelvinVectorSize + 2, 0)
            .noalias() = dF_dsigma.transpose() / G;
    }

    // G_54
    jacobian(2 * KelvinVectorSize + 2, 2 * KelvinVectorSize + 1) =
        -_mp.kappa(t, x)[0] * _mp.hardening_coefficient(t, x)[0] / G;

    // G_52, G_53, G_55 are zero
}

template <int DisplacementDim>
void SolidEhlers<DisplacementDim>::updateDamage(
    double const t, ProcessLib::SpatialPosition const& x,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
        material_state_variables)
{
    assert(dynamic_cast<MaterialStateVariables*>(&material_state_variables) !=
           nullptr);
    auto& state =
        static_cast<MaterialStateVariables&>(material_state_variables);

    // Default case of the rate problem. Updated below if volumetric plastic
    // strain rate is positive (dilatancy).
    state.damage.kappa_d = state.damage_prev.kappa_d;

    // Compute damage current step
    double const del_eps_p_V = state.eps_p_V - state.eps_p_V_prev;
    if (del_eps_p_V > 0)
    {
        double const h_d = _damage_properties->h_d(t, x)[0];
        double const del_eps_p_eff = state.eps_p_eff - state.eps_p_eff_prev;
        double const r_s = del_eps_p_eff / del_eps_p_V;

        // Brittleness decrease with confinement for the nonlinear flow rule.
        // ATTENTION: For linear flow rule -> constant brittleness.
        double x_s = 0;
        if (r_s < 1)
        {
            x_s = 1 + h_d * r_s * r_s;
        }
        else
        {
            x_s = 1 - 3 * h_d + 4 * h_d * std::sqrt(r_s);
        }
        state.damage.kappa_d += del_eps_p_eff / x_s;
    }

    double const alpha_d = _damage_properties->alpha_d(t, x)[0];
    double const beta_d = _damage_properties->beta_d(t, x)[0];

    // Update internal damage variable.
    state.damage.damage =
        (1 - beta_d) * (1 - std::exp(-state.damage.kappa_d / alpha_d));
}

/// Calculates the derivative of the residuals with respect to total
/// strain. Implementation fully implicit only.
template <int DisplacementDim>
ProcessLib::KelvinMatrixType<DisplacementDim> calculateDResidualDEps(
    double const K, double const G)
{
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;

    auto const& P_dev = Invariants::deviatoric_projection;
    auto const& P_sph = Invariants::spherical_projection;
    auto const& I = ProcessLib::KelvinMatrixType<DisplacementDim>::Identity();

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
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;
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

template <int DisplacementDim>
std::tuple<typename SolidEhlers<DisplacementDim>::KelvinVector,
           std::unique_ptr<
               typename MechanicsBase<DisplacementDim>::MaterialStateVariables>,
           typename SolidEhlers<DisplacementDim>::KelvinMatrix>
SolidEhlers<DisplacementDim>::integrateStress(
    double const t,
    ProcessLib::SpatialPosition const& x,
    double const dt,
    KelvinVector const& eps_prev,
    KelvinVector const& eps,
    KelvinVector const& sigma_prev,
    KelvinVector const& /*sigma_final*/,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
        material_state_variables)
{
    assert(dynamic_cast<MaterialStateVariables const*>(
               &material_state_variables) != nullptr);
    auto& _state =
        static_cast<MaterialStateVariables const&>(material_state_variables);

    MaterialStateVariables state(_state);
    state.setInitialConditions();

    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;

    // volumetric strain
    double const eps_V = Invariants::trace(eps);

    auto const& P_dev = Invariants::deviatoric_projection;
    // deviatoric strain
    KelvinVector const eps_D = P_dev * eps;

    // dimensionless stress/hydrostatic pressure
    double const G = _mp.G(t, x)[0];
    double const K = _mp.K(t, x)[0];

    KelvinVector sigma_eff_prev = sigma_prev;  // In case without damage the
                                               // effective value is same as the
                                               // previous one.
    if (_damage_properties)
    {
        // Compute sigma_eff from damage total stress sigma, which is given by
        // sigma_eff=sigma_prev / (1-damage)
        sigma_eff_prev = sigma_prev / (1 - state.damage_prev.damage);
    }
    KelvinVector sigma = predict_sigma<DisplacementDim>(G, K, sigma_eff_prev,
                                                        eps, eps_prev, eps_V);

    KelvinMatrix tangentStiffness;

    // update parameter
    double const k =
        calculateIsotropicHardening(_mp.kappa(t, x)[0],
                                    _mp.hardening_coefficient(t, x)[0],
                                    state.eps_p_eff);

    PhysicalStressWithInvariants<DisplacementDim> s{G * sigma};
    // Quit early if sigma is zero (nothing to do) or if we are still in elastic
    // zone.
    if (sigma.squaredNorm() == 0 ||
        yieldFunction<DisplacementDim>(s, _mp, t, x, k) < 0)
    {
        tangentStiffness.setZero();
        tangentStiffness.template topLeftCorner<3, 3>().setConstant(
            K - 2. / 3 * G);
        tangentStiffness.noalias() += 2 * G * KelvinMatrix::Identity();
    }
    else
    {
        JacobianMatrix jacobian;

        // Linear solver for the newton loop is required after the loop with the
        // same matrix. This saves one decomposition.
        Eigen::FullPivLU<JacobianMatrix> linear_solver;

        {  // Newton loop for return mapping calculation.
            auto const update_residual = [&](ResidualVectorType& residual) {

                KelvinVector const eps_p_D_dot =
                    (state.eps_p_D - state.eps_p_D_prev) / dt;
                double const eps_p_V_dot =
                    (state.eps_p_V - state.eps_p_V_prev) / dt;
                double const eps_p_eff_dot =
                    (state.eps_p_eff - state.eps_p_eff_prev) / dt;

                double const k = calculateIsotropicHardening(
                    _mp.kappa(t, x)[0],
                    _mp.hardening_coefficient(t, x)[0],
                    state.eps_p_eff);

                calculatePlasticResidual<DisplacementDim>(
                    t, x, eps_D, eps_V, s, state.eps_p_D, eps_p_D_dot,
                    state.eps_p_V, eps_p_V_dot, eps_p_eff_dot, state.lambda, k,
                    _mp, residual);
            };

            auto const update_jacobian = [&](JacobianMatrix& jacobian) {
                calculatePlasticJacobian<DisplacementDim>(dt, t, x, jacobian, s,
                                                          state.lambda, _mp);
            };

            auto const update_solution = [&](
                ResidualVectorType const& increment) {
                sigma.noalias() += increment.template segment<KelvinVectorSize>(
                    KelvinVectorSize * 0);
                s = PhysicalStressWithInvariants<DisplacementDim>{G * sigma};
                state.eps_p_D.noalias() +=
                    increment.template segment<KelvinVectorSize>(
                        KelvinVectorSize * 1);
                state.eps_p_V += increment(KelvinVectorSize * 2);
                state.eps_p_eff += increment(KelvinVectorSize * 2 + 1);
                state.lambda += increment(KelvinVectorSize * 2 + 2);
            };

            auto newton_solver = NumLib::NewtonRaphson<
                decltype(linear_solver), JacobianMatrix,
                decltype(update_jacobian), ResidualVectorType,
                decltype(update_residual), decltype(update_solution)>(
                linear_solver, update_jacobian, update_residual,
                update_solution, _nonlinear_solver_parameters);

            auto const success_iterations = newton_solver.solve(jacobian);

            if (!success_iterations)
                return {sigma, nullptr, tangentStiffness};

            // If the Newton loop didn't run, the linear solver will not be
            // initialized.
            // This happens usually for the first iteration of the first
            // timestep.
            if (*success_iterations == 0)
                linear_solver.compute(jacobian);
        }

        // Calculate residual derivative w.r.t. strain
        Eigen::Matrix<double, JacobianResidualSize, KelvinVectorSize,
                      Eigen::RowMajor>
            dresidual_deps =
                Eigen::Matrix<double, JacobianResidualSize, KelvinVectorSize,
                              Eigen::RowMajor>::Zero();
        dresidual_deps.template block<KelvinVectorSize, KelvinVectorSize>(0, 0)
            .noalias() = calculateDResidualDEps<DisplacementDim>(K, G);

        if (_damage_properties)
            updateDamage(t, x, state);

        // Extract consistent tangent.
        tangentStiffness =
            _mp.G(t, x)[0] *
            linear_solver.solve(-dresidual_deps)
                .template block<KelvinVectorSize, KelvinVectorSize>(0, 0);
    }

    if (_damage_properties)
        return {_mp.G(t, x)[0] * sigma * (1 - state.damage.damage),
                std::unique_ptr<typename MechanicsBase<
                    DisplacementDim>::MaterialStateVariables>{
                    new MaterialStateVariables{state}},
                std::move(tangentStiffness)};

    return {
        _mp.G(t, x)[0] * sigma,
        std::unique_ptr<
            typename MechanicsBase<DisplacementDim>::MaterialStateVariables>{
            new MaterialStateVariables{state}},
        tangentStiffness};
}

}  // namespace Ehlers
}  // namespace Solids
}  // namespace MaterialLib
