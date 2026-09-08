// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>
#include <string>

namespace MaterialPropertyLib
{
/// State the drift-flux closure is evaluated at. The closure, its quadratic
/// coefficients, its residual, the diagnostics of a state it cannot solve, and
/// the slip parameter of the resulting mixture all take the same state, so a
/// caller cannot describe one state and solve another.
struct DriftFluxState
{
    /// The vapour mass fraction, dimensionless.
    double dryness;
    /// In kg/m^3.
    double vapour_water_density;
    /// In kg/m^3.
    double liquid_water_density;
    /// The mixture velocity in m/s.
    double v_mix;
    /// The profile parameter, dimensionless, see
    /// driftFluxProfileParameter().
    double C_0;
    /// The drift flux velocity in m/s, see driftFluxVelocity().
    double u_gu;
};

/// Coefficients of the residual divided by the liquid density,
/// \f$R / \rho_l = a \alpha^2 + b \alpha + c\f$. All three are mass fluxes,
/// in kg/(m^2 s), since the void fraction is dimensionless.
struct VoidFractionQuadratic
{
    double a;
    double b;
    double c;
};

/// Profile parameter \f$C_0\f$ of the Rouhani-Axelsson drift-flux closure,
/// the flow-weighted ratio of the cross-sectional averages that accounts for
/// the non-uniform void and velocity profiles over the well cross-section.
/// \param dryness  the vapour mass fraction, dimensionless.
/// \return The profile parameter, dimensionless.
double driftFluxProfileParameter(double const dryness);

/// Coefficients of the drift-flux closure written as a quadratic equation in
/// the void fraction, see computeVapourVoidFraction().
///
/// \param state  the state of the closure.
/// \return The coefficients, each in kg/(m^2 s).
VoidFractionQuadratic voidFractionQuadratic(DriftFluxState const& state);

/// Residual of the drift-flux closure. Kept for testing and for the residuals
/// reported with a non-solvable state; the void fraction itself is computed in
/// closed form.
///
/// \param alpha  the vapour void fraction, dimensionless.
/// \param state  the state of the closure.
/// \return The residual in kg^2/(m^5 s), the mass flux of the closure times
/// the liquid density it is written with.
double voidFractionResidual(double const alpha, DriftFluxState const& state);

/// State of the drift-flux closure and its residuals at the ends of the
/// admissible interval \f$[0, x / S]\f$, formatted for the message of an
/// assembly aborted because computeVapourVoidFraction() found no admissible
/// void fraction. The residuals are reported for the closure and the interval
/// that computeVapourVoidFraction() actually solves on, that is with the drift
/// aligned with the mixture flow and up to \f$x / S\f$ rather than one.
///
/// The caller prepends its own context, the process or boundary condition, the
/// element, and the primary variables it has at hand.
///
/// \param state  the state of the closure, reported with its residuals.
std::string voidFractionClosureDiagnostics(DriftFluxState const& state);

/// Drift flux velocity \f$u_{gu}\f$ of the Rouhani-Axelsson closure, from the
/// surface tension of the IAPWS correlation, see Cooper, J. R., and R. B.
/// Dooley. "IAPWS release on surface tension of ordinary water substance."
/// International Association for the Properties of Water and Steam (1994).
///
/// Both the surface tension correlation and the buoyancy driving the drift are
/// capped at zero at the critical point, where the two phases become
/// identical: beyond it the reduced temperature is negative and the phase
/// densities cross, and std::pow of a negative base with the non-integer
/// exponents of the two correlations is NaN, which would travel into the void
/// fraction closure. Newton iterates do reach that range, so both are capped
/// rather than assumed positive.
///
/// \param dryness                the vapour mass fraction, dimensionless.
/// \param temperature            in K.
/// \param vapour_water_density   in kg/m^3.
/// \param liquid_water_density   in kg/m^3.
/// \return The drift flux velocity in m/s.
double driftFluxVelocity(double const dryness, double const temperature,
                         double const vapour_water_density,
                         double const liquid_water_density);

/// Drift flux velocity aligned with the mixture flow,
/// \f$\operatorname{sign}(v)\, u_{gu}\f$. The Rouhani-Axelsson closure is
/// derived for co-current flow and has no admissible solution for backflow
/// with a drift velocity that is fixed in the gravity frame, see
/// computeVapourVoidFraction() for the reasoning.
///
/// Every use of the drift flux velocity next to a void fraction obtained from
/// that closure, in particular the slip momentum term of the mixture, has to
/// use the same aligned value: the void fraction is even in the mixture
/// velocity, so a raw drift flux velocity leaves the slip term inconsistent
/// with it and lets the term vanish at the finite backflow velocity
/// \f$v = -u_{gu} / (C_0 - 1)\f$.
double alignedDriftFluxVelocity(double const u_gu, double const v_mix);

/// State of the drift-flux closure for a two-phase water mixture, assembled
/// from the quantities a local assembler has at an integration point. The
/// profile parameter and the drift flux velocity are not independent inputs:
/// they follow from the same dryness, temperature and phase densities through
/// driftFluxProfileParameter() and driftFluxVelocity(), and the drift has to
/// be the one aligned with the mixture flow, see alignedDriftFluxVelocity().
/// Composing them here keeps every caller on the same closure.
///
/// \param dryness                the vapour mass fraction, dimensionless.
/// \param temperature            in K.
/// \param vapour_water_density   in kg/m^3.
/// \param liquid_water_density   in kg/m^3.
/// \param v_mix                  the mixture velocity in m/s.
/// \return The state, with its profile parameter and its aligned drift flux
/// velocity.
DriftFluxState driftFluxState(double const dryness, double const temperature,
                              double const vapour_water_density,
                              double const liquid_water_density,
                              double const v_mix);

/// Vapour void fraction of the Rouhani-Axelsson drift-flux closure, see
/// Rouhani, Z., and E. Axelsson. "Calculation of volume void fraction in a
/// subcooled and quality region." International Journal of Heat and Mass
/// Transfer 17 (1970): 383-393.
///
/// With the mixture density \f$\rho_m = \alpha \rho_v + (1 - \alpha) \rho_l\f$
/// and the mass flux \f$G = \rho_m v\f$ the closure reads
/// \f[
///     G (x - \alpha S) - \alpha \rho_v u_{gu} = 0, \quad
///     S = C_0 \left(x + (1 - x) \frac{\rho_v}{\rho_l}\right),
/// \f]
/// which is a quadratic equation in \f$\alpha\f$ because \f$\rho_m\f$ depends
/// on \f$\alpha\f$ linearly.
///
/// For co-current upflow, \f$v > 0\f$ and \f$0 < x < 1\f$, the residual is
/// positive at \f$\alpha = 0\f$ and negative at \f$\alpha = 1\f$, hence a root
/// exists, and it is the only root in the admissible interval \f$[0, \min(1, x
/// / S)]\f$. The bound \f$x / S\f$ is the homogeneous void fraction divided by
/// the profile parameter and is at most one as long as \f$C_0 \ge 1\f$, which
/// is what driftFluxProfileParameter() gives; the minimum keeps the returned
/// volume fraction below one for a profile parameter below one as well.
/// As long as the phase densities are ordered, \f$\rho_l > \rho_v\f$, the
/// parabola opens upwards and the second root is larger than one. Where the
/// two densities cross, in the immediate vicinity of the critical point, it
/// opens downwards and both roots can lie below one, which is why the root is
/// selected by testing the admissible interval rather than by its position.
///
/// The correlation is derived for co-current flow. Written with a drift
/// velocity that is fixed in the gravity frame it has no admissible solution
/// for backflow, \f$v < 0\f$: the denominator \f$S G + \rho_v u_{gu}\f$ of the
/// closure changes sign, so the void fraction has a pole and leaves \f$[0,
/// 1]\f$. The drift velocity is therefore aligned with the mixture flow,
/// \f$u_{gu} \to \operatorname{sign}(v)\, u_{gu}\f$, which makes the residual
/// odd under \f$v \to -v\f$ and hence the void fraction an even function of
/// the mixture velocity. Existence and uniqueness above then hold for every
/// mixture velocity, and \f$\alpha \sim |v|\f$ near \f$v = 0\f$, that is the
/// void fraction is continuous but not differentiable where the flow reverses.
/// They hold for the exact closure; the arithmetic that implements it reports
/// no admissible root once the mixture velocity is small enough that the
/// coefficients of the quadratic underflow, below some \f$10^{-160}\f$ m/s,
/// which is answered by a retreat like any other state without a root. A
/// mixture at rest, \f$v = 0\f$ exactly, is resolved by continuity instead,
/// see below.
///
/// The kink at \f$v = 0\f$ could be removed by scaling the drift velocity with
/// \f$v / \sqrt{v^2 + u_{gu}^2}\f$, using the drift velocity itself as the
/// velocity scale so that no tuning parameter is introduced. This is not done
/// because it changes the void fraction in the whole band \f$|v| \lesssim
/// u_{gu}\f$; what the kink does to the convergence of the global Newton
/// solver has not been investigated.
///
/// Solves the closure in closed form and returns the root in the admissible
/// interval, or no value if the closure has no admissible solution for the
/// given state.
///
/// A mixture at rest without drift, \f$v = 0\f$ and \f$u_{gu} = 0\f$, leaves
/// every coefficient of the quadratic zero and the closure is satisfied by
/// every void fraction. It is resolved by continuity: without drift the
/// closure reads \f$x = \alpha S\f$ for every non-zero mixture velocity,
/// independently of it, so the limit is the profile slip alone, \f$\alpha = x
/// / S\f$.
///
/// Single phase states, \f$x \le 0\f$ and \f$x \ge 1\f$, return \f$0\f$ and
/// \f$1\f$ without looking at the state, in particular without looking at the
/// density of the absent phase.
///
/// For a two-phase state it throws NumLib::AssemblyException for a
/// non-positive phase density and for a non-finite dryness, mixture velocity,
/// drift flux velocity, or profile parameter: all of those follow the solution
/// iterate, so the assembly is aborted rather than the run ended.
///
/// What the abort leads to is up to the nonlinear solver, and both of them
/// answer it the same way: the Newton and the Picard solver each catch
/// NumLib::AssemblyException, end the nonlinear iteration and let the time
/// stepping repeat the step with a smaller time step size, with the state of
/// the closure reported.
///
/// \param state  the state of the closure.
/// \return The vapour void fraction, dimensionless, or no value if the closure
/// has no admissible solution for that state.
std::optional<double> computeVapourVoidFraction(DriftFluxState const& state);

/// Slip parameter \f$\gamma\f$ of the two-phase mixture, the momentum flux
/// carried by the relative motion of the phases, see Akbar, Somaieh, N.
/// Fathianpour, and Rafid Al Khoury. "A finite element model for high enthalpy
/// two-phase flow in geothermal wellbores." Renewable Energy 94 (2016):
/// 223-236.
/// \f[
///     \gamma = \frac{\alpha \rho_l \rho_v \rho_m}{(1 - \alpha)
///              \left(\alpha C_0 \rho_v + (1 - \alpha C_0) \rho_l\right)^2}
///              \left((C_0 - 1) v + u_{gu}\right)^2
/// \f]
/// with the mixture density \f$\rho_m = \alpha \rho_v + (1 - \alpha)
/// \rho_l\f$.
///
/// A void fraction of one, which computeVapourVoidFraction() returns at a
/// dryness of one, would divide by zero. There is no liquid phase left to slip
/// against and the limit is zero slip: both the liquid fraction \f$1 -
/// \alpha\f$ and the bracket \f$(C_0 - 1) v + u_{gu}\f$ vanish linearly in
/// \f$1 - x\f$, and the bracket enters squared, so the quotient is of the
/// order of \f$1 - x\f$.
///
/// The dryness of \c state does not enter; the void fraction it leads to is
/// passed instead.
///
/// \param alpha  the vapour void fraction, dimensionless.
/// \param state  the state of the closure the void fraction was solved from.
/// \return The slip parameter in Pa; it is a momentum flux and enters the
/// momentum balance next to \f$\rho_m v^2\f$.
double mixtureSlipParameter(double const alpha, DriftFluxState const& state);
}  // namespace MaterialPropertyLib
