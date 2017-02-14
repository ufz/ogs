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
#include "NumLib/NewtonRaphson.h"

namespace MaterialLib
{
namespace Solids
{
namespace Ehlers
{
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
double yieldFunction(
    PhysicalStressWithInvariants<DisplacementDim> const& s,
    typename SolidEhlers<DisplacementDim>::MaterialProperties const& mp,
    double const t, ProcessLib::SpatialPosition const& x)
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
           beta * s.I_1 + epsilon * I_1_squared - mp.k;
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
    typename SolidEhlers<DisplacementDim>::MaterialProperties const& _mp,
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
        yieldFunction<DisplacementDim>(s, _mp, t, x) / G;
}

/// Special product of \p v with itself: \f$v \odot v\f$.
/// v is given in \p Kelvin mapping.
template <int DisplacementDim>
ProcessLib::KelvinMatrixType<DisplacementDim> s_odot_s(
    ProcessLib::KelvinVectorType<DisplacementDim> const& v)
{
    ProcessLib::KelvinMatrixType<DisplacementDim> result;

    result(0, 0) = v(0) * v(0);
    result(0, 1) = result(1, 0) = v(3) * v(3) / 2.;
    result(0, 2) = result(2, 0) = v(5) * v(5) / 2.;
    result(0, 3) = result(3, 0) = v(0) * v(3);
    result(0, 4) = result(4, 0) = v(3) * v(5) / std::sqrt(2.);
    result(0, 5) = result(5, 0) = v(0) * v(5);

    result(1, 1) = v(1) * v(1);
    result(1, 2) = result(2, 1) = v(4) * v(4) / 2.;
    result(1, 3) = result(3, 1) = v(3) * v(1);
    result(1, 4) = result(4, 1) = v(1) * v(4);
    result(1, 5) = result(5, 1) = v(3) * v(4) / std::sqrt(2.);

    result(2, 2) = v(2) * v(2);
    result(2, 3) = result(3, 2) = v(5) * v(4) / std::sqrt(2.);
    result(2, 4) = result(4, 2) = v(4) * v(2);
    result(2, 5) = result(5, 2) = v(5) * v(2);

    result(3, 3) = v(0) * v(1) + v(3) * v(3) / 2.;
    result(3, 4) = result(4, 3) =
        v(3) * v(4) / 2. + v(5) * v(1) / std::sqrt(2.);
    result(3, 5) = result(5, 3) =
        v(0) * v(4) / std::sqrt(2.) + v(3) * v(5) / 2.;

    result(4, 4) = v(1) * v(2) + v(4) * v(4) / 2.;
    result(4, 5) = result(5, 4) =
        v(3) * v(2) / std::sqrt(2.) + v(5) * v(4) / 2.;

    result(5, 5) = v(0) * v(2) + v(5) * v(5) / 2.;
    return result;
}

template <int DisplacementDim>
void calculatePlasticJacobian(
    double const dt,
    double const t,
    ProcessLib::SpatialPosition const& x,
    typename SolidEhlers<DisplacementDim>::JacobianMatrix& jacobian,
    PhysicalStressWithInvariants<DisplacementDim> const& s,
    double const lambda,
    typename SolidEhlers<DisplacementDim>::MaterialProperties const& _mp)
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
        theta * P_dev * s_odot_s<DisplacementDim>(sigma_D_inverse) * P_dev +
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

template <int DisplacementDim>
void SolidEhlers<DisplacementDim>::MaterialProperties::
    calculateIsotropicHardening(double const t,
                                ProcessLib::SpatialPosition const& x,
                                double const eps_p_eff)
{
    k = kappa(t, x)[0] * (1. + eps_p_eff * hardening_coefficient(t, x)[0]);
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
bool SolidEhlers<DisplacementDim>::computeConstitutiveRelation(
    double const t,
    ProcessLib::SpatialPosition const& x,
    double const dt,
    KelvinVector const& eps_prev,
    KelvinVector const& eps,
    KelvinVector const& sigma_prev,
    KelvinVector& sigma_final,
    KelvinMatrix& C,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
        material_state_variables)
{
    assert(dynamic_cast<MaterialStateVariables*>(&material_state_variables) !=
           nullptr);
    MaterialStateVariables& _state =
        static_cast<MaterialStateVariables&>(material_state_variables);
    _state.setInitialConditions();

    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;
    C.setZero();

    // volumetric strain
    double const eps_V = Invariants::trace(eps);

    auto const& P_dev = Invariants::deviatoric_projection;
    // deviatoric strain
    KelvinVector const eps_D = P_dev * eps;

    // dimensionless stress/hydrostatic pressure
    double const G = _mp.G(t, x)[0];
    double const K = _mp.K(t, x)[0];

    KelvinVector sigma =
        predict_sigma<DisplacementDim>(G, K, sigma_prev, eps, eps_prev, eps_V);

    // update parameter
    _mp.calculateIsotropicHardening(t, x, _state.eps_p_eff);

    PhysicalStressWithInvariants<DisplacementDim> s{G * sigma};
    // Quit early if sigma is zero (nothing to do) or if we are still in elastic
    // zone.
    if (sigma.squaredNorm() == 0 ||
        yieldFunction<DisplacementDim>(s, _mp, t, x) < 0)
    {
        C.setZero();
        C.template topLeftCorner<3, 3>().setConstant(K - 2. / 3 * G);
        C.noalias() += 2 * G * KelvinMatrix::Identity();
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
                    (_state.eps_p_D - _state.eps_p_D_prev) / dt;
                double const eps_p_V_dot =
                    (_state.eps_p_V - _state.eps_p_V_prev) / dt;
                double const eps_p_eff_dot =
                    (_state.eps_p_eff - _state.eps_p_eff_prev) / dt;
                calculatePlasticResidual<DisplacementDim>(
                    t, x, eps_D, eps_V, s, _state.eps_p_D, eps_p_D_dot,
                    _state.eps_p_V, eps_p_V_dot, eps_p_eff_dot, _state.lambda,
                    _mp, residual);
            };

            auto const update_jacobian = [&](JacobianMatrix& jacobian) {
                calculatePlasticJacobian<DisplacementDim>(dt, t, x, jacobian, s,
                                                          _state.lambda, _mp);
            };

            auto const update_solution = [&](
                ResidualVectorType const& increment) {
                sigma.noalias() += increment.template segment<KelvinVectorSize>(
                    KelvinVectorSize * 0);
                s = PhysicalStressWithInvariants<DisplacementDim>{G * sigma};
                _state.eps_p_D.noalias() +=
                    increment.template segment<KelvinVectorSize>(
                        KelvinVectorSize * 1);
                _state.eps_p_V += increment(KelvinVectorSize * 2);
                _state.eps_p_eff += increment(KelvinVectorSize * 2 + 1);
                _state.lambda += increment(KelvinVectorSize * 2 + 2);

                _mp.calculateIsotropicHardening(t, x, _state.eps_p_eff);
            };

            // TODO Make the following choice of maximum iterations and
            // convergence criteria available from the input file configuration:
            int const maximum_iterations(100);
            double const tolerance(1e-14);

            auto newton_solver = NumLib::NewtonRaphson<
                decltype(linear_solver), JacobianMatrix,
                decltype(update_jacobian), ResidualVectorType,
                decltype(update_residual), decltype(update_solution)>(
                linear_solver, update_jacobian, update_residual,
                update_solution, maximum_iterations, tolerance);

            auto const success_iterations = newton_solver.solve(jacobian);

            if (!success_iterations)
                return false;

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

        // Extract consistent tangent.
        C.noalias() =
            _mp.G(t, x)[0] *
            linear_solver.solve(-dresidual_deps)
                .template block<KelvinVectorSize, KelvinVectorSize>(0, 0);
    }

    // Update sigma.
    sigma_final.noalias() = G * sigma;
    return true;
}

}  // namespace Ehlers
}  // namespace Solids
}  // namespace MaterialLib
