/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Adsorption.h"

#include "BaseLib/Logging.h"
#include "MaterialLib/PhysicalConstant.h"

namespace
{
// Evaluate adsorbtion potential A
double getPotential(const double p_Ads, const double T_Ads, const double M_Ads)
{
    return MaterialLib::PhysicalConstant::IdealGasConstant * T_Ads *
           std::log(
               Adsorption::AdsorptionReaction::getEquilibriumVapourPressure(
                   T_Ads) /
               p_Ads) /
           (M_Ads * 1.e3);  // in kJ/kg = J/g
}

constexpr double k_rate = 6.0e-3;  // to be specified

template <typename T>
T square(const T& v)
{
    return v * v;
}
}  // namespace

namespace Adsorption
{
// Saturation pressure for water used in Nunez
double AdsorptionReaction::getEquilibriumVapourPressure(const double T_Ads)
{
    // critical T and p
    const double Tc = 647.3;    // K
    const double pc = 221.2e5;  // Pa
    // dimensionless T
    const double Tr = T_Ads / Tc;
    const double theta = 1. - Tr;
    // empirical constants
    const double c[] = {-7.69123, -26.08023, -168.17065, 64.23285, -118.96462,
                        4.16717,  20.97506,  1.0e9,      6.0};
    const double K[] = {
        c[0] * theta + c[1] * std::pow(theta, 2) + c[2] * std::pow(theta, 3) +
            c[3] * std::pow(theta, 4) + c[4] * std::pow(theta, 5),
        1. + c[5] * theta + c[6] * std::pow(theta, 2)};

    const double exponent =
        K[0] / (K[1] * Tr) - theta / (c[7] * std::pow(theta, 2) + c[8]);
    return pc * std::exp(exponent);  // in Pa
}

// Evaporation enthalpy of water from Nunez
double AdsorptionReaction::getEvaporationEnthalpy(double T_Ads)  // in kJ/kg
{
    T_Ads -= 273.15;
    if (T_Ads <= 10.)
    {
        const double c[] = {2.50052e3,   -2.1068,     -3.57500e-1,
                            1.905843e-1, -5.11041e-2, 7.52511e-3,
                            -6.14313e-4, 2.59674e-5,  -4.421e-7};
        double hv = 0.;
        for (size_t i = 0; i < sizeof(c) / sizeof(c[0]); i++)
        {
            hv += c[i] * std::pow(T_Ads, i);
        }
        return hv;
    }
    if (T_Ads <= 300.)
    {
        const double c[] = {2.50043e3,  -2.35209,    1.91685e-4,  -1.94824e-5,
                            2.89539e-7, -3.51199e-9, 2.06926e-11, -6.4067e-14,
                            8.518e-17,  1.558e-20,   -1.122e-22};
        double hv = 0.;
        for (size_t i = 0; i < sizeof(c) / sizeof(c[0]); i++)
        {
            hv += c[i] * std::pow(T_Ads, i);
        }
        return hv;
    }
    const double c[] = {2.99866e3, -3.1837e-3,  -1.566964e1,
                        -2.514e-6, 2.045933e-2, 1.0389e-8};
    return ((c[0] + c[2] * T_Ads + c[4] * std::pow(T_Ads, 2)) /
            (1. + c[1] * T_Ads + c[3] * std::pow(T_Ads, 2) +
             c[5] * std::pow(T_Ads, 3)));
}

double AdsorptionReaction::getMolarFraction(double xm, double M_this,
                                            double M_other)
{
    return M_other * xm / (M_other * xm + M_this * (1.0 - xm));
}

double AdsorptionReaction::dMolarFraction(double xm, double M_this,
                                          double M_other)
{
    return M_other * M_this / square(M_other * xm + M_this * (1.0 - xm));
}

double AdsorptionReaction::getReactionRate(const double p_Ads,
                                           const double T_Ads,
                                           const double M_Ads,
                                           const double loading) const
{
    const double A = getPotential(p_Ads, T_Ads, M_Ads);
    const double C_eq =
        std::max(0., getAdsorbateDensity(T_Ads) * characteristicCurve(A));

    return k_rate * (C_eq - loading);  // scaled with mass fraction
                                       // this the rate in terms of loading!
}

double AdsorptionReaction::getLoading(const double rho_curr,
                                      const double rho_dry)
{
    return rho_curr / rho_dry - 1.0;
}

// Calculate sorption entropy
double AdsorptionReaction::getEntropy(const double T_Ads, const double A) const
{
    const double epsilon = 1.0e-8;

    //* // This change will also change simulation results.
    const double W_p = characteristicCurve(A + epsilon);
    const double W_m = characteristicCurve(A - epsilon);
    const double dAdlnW = 2.0 * epsilon / (std::log(W_p / W_m));
    // */

    if (W_p <= 0.0 || W_m <= 0.0)
    {
        ERR("characteristic curve in negative region (W-, W+): {:g}, {:g}", W_m,
            W_p);
        return 0.0;
    }

    return dAdlnW * getAlphaT(T_Ads);
}

// Calculate sorption enthalpy
double AdsorptionReaction::getEnthalpy(const double p_Ads, const double T_Ads,
                                       const double M_Ads) const
{
    // TODO [CL] consider using A as obtained from current loading (needs
    // inverse CC A(W)) instead of p_Vapour, T_Vapour
    const double A = getPotential(p_Ads, T_Ads, M_Ads);

    return (getEvaporationEnthalpy(T_Ads) + A - T_Ads * getEntropy(T_Ads, A)) *
           1000.0;  // in J/kg
}

double AdsorptionReaction::getEquilibriumLoading(const double p_Ads,
                                                 const double T_Ads,
                                                 const double M_Ads) const
{
    const double A = getPotential(p_Ads, T_Ads, M_Ads);
    return getAdsorbateDensity(T_Ads) * characteristicCurve(A);
}

}  // namespace Adsorption
