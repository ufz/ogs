// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace MaterialPropertyLib
{
/// Profile parameter \f$C_0\f$ of the Rouhani-Axelsson drift-flux closure,
/// the flow-weighted ratio of the cross-sectional averages that accounts for
/// the non-uniform void and velocity profiles over the well cross-section.
/// \param dryness  the vapour mass fraction, dimensionless.
double driftFluxProfileParameter(double const dryness);

/// Drift flux velocity \f$u_{gu}\f$ of the Rouhani-Axelsson closure, from the
/// surface tension of the IAPWS correlation, see Cooper, J. R., and R. B.
/// Dooley. "IAPWS release on surface tension of ordinary water substance."
/// International Association for the Properties of Water and Steam (1994).
double driftFluxVelocity(double const dryness, double const temperature,
                         double const vapour_water_density,
                         double const liquid_water_density);

/// Vapour void fraction of the Rouhani-Axelsson drift-flux closure, see
/// Rouhani, Z., and E. Axelsson. "Calculation of volume void fraction in a
/// subcooled and quality region." International Journal of Heat and Mass
/// Transfer 17 (1970): 383-393.
///
/// Solved with a local Newton iteration started from zero. Non-convergence is
/// warned about and the last iterate is returned.
double computeVapourVoidFraction(double const dryness,
                                 double const vapour_water_density,
                                 double const liquid_water_density,
                                 double const v_mix, double const C_0,
                                 double const u_gu);

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
/// \param alpha                  the vapour void fraction, dimensionless.
/// \param vapour_water_density   in kg/m^3.
/// \param liquid_water_density   in kg/m^3.
/// \param v_mix                  the mixture velocity in m/s.
/// \param C_0                    the profile parameter, dimensionless.
/// \param u_gu                   the drift flux velocity in m/s.
double mixtureSlipParameter(double const alpha,
                            double const vapour_water_density,
                            double const liquid_water_density,
                            double const v_mix, double const C_0,
                            double const u_gu);
}  // namespace MaterialPropertyLib
