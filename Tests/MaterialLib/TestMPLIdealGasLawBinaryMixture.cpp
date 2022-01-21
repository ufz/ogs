/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include "MaterialLib/MPL/Properties/IdealGasLawBinaryMixture.h"
#include "MaterialLib/PhysicalConstant.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

double molarFractionFantasy(const double pGR, const double pGR_min,
                            const double pGR_max, const double pCap,
                            const double pCap_min, const double pCap_max,
                            const double T, const double T_min,
                            const double T_max)
{
    auto const dpGR = pGR_max - pGR_min;
    auto const dpCap = pCap_max - pCap_min;
    auto const dT = T_max - T_min;

    auto const a = 1. / (6. * std::pow(dpGR, 2.));
    auto const b = 1. / (6. * dpGR);
    auto const c = 1. / (6. * std::pow(dpCap, 2.));
    auto const d = 1. / (6. * (dpCap));
    auto const e = 1. / (6. * std::pow(dT, 2.));
    auto const f = 1. / (6. * (dT));

    return a * (pGR - pGR_min) * (pGR - pGR_min) + b * (pGR - pGR_min) +
           c * (pCap - pCap_min) * (pCap - pCap_min) + d * (pCap - pCap_min) +
           e * (T - T_min) * (T - T_min) + f * (T - T_min);
}

std::array<double, 3> molarMassFantasyDerivatives(
    const double pGR, const double pGR_min, const double pGR_max,
    const double pCap, const double pCap_min, const double pCap_max,
    const double T, const double T_min, const double T_max)
{
    auto const dpGR = pGR_max - pGR_min;
    auto const dpCap = pCap_max - pCap_min;
    auto const dT = T_max - T_min;

    auto const a = 1. / (6. * std::pow(dpGR, 2.));
    auto const b = 1. / (6. * dpGR);
    auto const c = 1. / (6. * std::pow(dpCap, 2.));
    auto const d = 1. / (6. * (dpCap));
    auto const e = 1. / (6. * std::pow(dT, 2.));
    auto const f = 1. / (6. * (dT));

    const double d_dpGR = 2 * a * (pGR - pGR_min) + b;
    const double d_dpCap = 2 * c * (pCap - pCap_min) + d;
    const double d_dT = 2 * e * (T - T_min) + f;

    return {d_dpGR, d_dpCap, d_dT};
}

TEST(MaterialPropertyLib, IdealGasLawBinaryMixture)
{
    const double MC = 0.028949;
    const double MW = 0.018052;

    auto const density_model = MPL::IdealGasLawBinaryMixture("density");

    enum
    {
        gas_pressure,
        capillary_pressure,
        temperature
    };
    enum
    {
        gas_phase_density,
        vapour_density,
        air_density
    };

    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;

    vars[static_cast<int>(MPL::Variable::molar_mass_vapour)] = MW;
    vars[static_cast<int>(MPL::Variable::molar_mass)] = MC;

    // Density derivatives of gases are dependent on pressure and temperature.
    // In the case of binary gases, there is also the dependence of the
    // composition (e.g. molar fractions).

    // The molar fraction of the vapour component in the gas phase is
    // represented here by a fantasy function. In Reality, this is actually a
    // physical quantity that depends on the temperature-dependent vapour
    // pressure, the capillary pressure-dependent correction factor for curved
    // menisci and the pressure of the gas phase.

    // In this test, however, the function is merely a polynomial,
    // which over the entire parameter space pGR-pCap-T assumes
    // values in the range [0,1].

    // Boundaries of the parameter space
    auto const pGR_min = 80000.;
    auto const pGR_max = 320000.;
    auto const pCap_min = 100.;
    auto const pCap_max = 1.e8;
    auto const T_min = 270.;
    auto const T_max = 520.;

    // travel through the entire parameter space
    for (double pGR = pGR_min; pGR <= pGR_max; pGR += 200.)
    {
        for (double pCap = pCap_min; pCap <= pCap_max; pCap *= 2.)
        {
            for (double T = T_min; T <= T_max; T += 0.2)
            {
                vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR;
                vars[static_cast<int>(MPL::Variable::temperature)] = T;

                // Molar fraction fantasy function
                auto const xnWG =
                    molarFractionFantasy(pGR, pGR_min, pGR_max, pCap, pCap_min,
                                         pCap_max, T, T_min, T_max);

                // ... and its derivatives
                auto const dxnWG = molarMassFantasyDerivatives(
                    pGR, pGR_min, pGR_max, pCap, pCap_min, pCap_max, T, T_min,
                    T_max);

                vars[static_cast<int>(MPL::Variable::molar_fraction)] = xnWG;
                vars[static_cast<int>(MPL::Variable::dxn_dpGR)] =
                    dxnWG[gas_pressure];
                vars[static_cast<int>(MPL::Variable::dxn_dpCap)] =
                    dxnWG[capillary_pressure];
                vars[static_cast<int>(MPL::Variable::dxn_dT)] =
                    dxnWG[temperature];

                auto const density = std::get<Eigen::Matrix<double, 3, 1>>(
                    density_model.value(vars, pos, t, dt));

                auto const rhoGR = density[gas_phase_density];
                auto const rhoCGR = density[air_density];
                auto const rhoWGR = density[vapour_density];

                const double R =
                    MaterialLib::PhysicalConstant::IdealGasConstant;

                auto const rhoCGR_ = (1. - xnWG) * pGR * MC / R / T;
                auto const rhoWGR_ = xnWG * pGR * MW / R / T;
                auto const rhoGR_ = rhoCGR_ + rhoWGR_;

                ASSERT_NEAR(rhoGR, rhoGR_, 1.e-18);
                ASSERT_NEAR(rhoCGR, rhoCGR_, 1.e-18);
                ASSERT_NEAR(rhoWGR, rhoWGR_, 1.e-18);

                auto const density_derivative =
                    std::get<Eigen::Matrix<double, 3, 3>>(density_model.dValue(
                        vars, MPL::Variable::phase_pressure, pos, t, dt));

                // Check derivatives via central differences
                // Pertubation of gas pressure
                auto const eps_pGR = 1.;

                auto const pGR_plus = pGR + eps_pGR;
                auto const xnWG_plus =
                    molarFractionFantasy(pGR_plus, pGR_min, pGR_max, pCap,
                                         pCap_min, pCap_max, T, T_min, T_max);

                vars[static_cast<int>(MPL::Variable::phase_pressure)] =
                    pGR_plus;
                vars[static_cast<int>(MPL::Variable::molar_fraction)] =
                    xnWG_plus;

                auto density_plus = std::get<Eigen::Matrix<double, 3, 1>>(
                    density_model.value(vars, pos, t, dt));

                auto const pGR_minus = pGR - eps_pGR;
                vars[static_cast<int>(MPL::Variable::phase_pressure)] =
                    pGR_minus;
                vars[static_cast<int>(MPL::Variable::molar_fraction)] =
                    molarFractionFantasy(pGR_minus, pGR_min, pGR_max, pCap,
                                         pCap_min, pCap_max, T, T_min, T_max);

                auto density_minus = std::get<Eigen::Matrix<double, 3, 1>>(
                    density_model.value(vars, pos, t, dt));

                for (std::size_t d = gas_phase_density; d <= air_density; d++)
                {
                    auto const drho_dpGR =
                        (density_plus[d] - density_minus[d]) / (2 * eps_pGR);

                    ASSERT_NEAR(drho_dpGR, density_derivative(d, gas_pressure),
                                1.e-10);
                }

                // Pertubation of capillary pressure
                auto const eps_pCap = 1.;
                auto const pCap_plus = pCap + eps_pCap;

                vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR;
                vars[static_cast<int>(MPL::Variable::capillary_pressure)] =
                    pCap_plus;
                vars[static_cast<int>(MPL::Variable::molar_fraction)] =
                    molarFractionFantasy(pGR, pGR_min, pGR_max, pCap_plus,
                                         pCap_min, pCap_max, T, T_min, T_max);

                density_plus = std::get<Eigen::Matrix<double, 3, 1>>(
                    density_model.value(vars, pos, t, dt));

                auto const pCap_minus = pCap - eps_pCap;
                vars[static_cast<int>(MPL::Variable::capillary_pressure)] =
                    pCap_minus;
                vars[static_cast<int>(MPL::Variable::molar_fraction)] =
                    molarFractionFantasy(pGR, pGR_min, pGR_max, pCap_minus,
                                         pCap_min, pCap_max, T, T_min, T_max);

                density_minus = std::get<Eigen::Matrix<double, 3, 1>>(
                    density_model.value(vars, pos, t, dt));

                for (std::size_t d = gas_phase_density; d <= air_density; d++)
                {
                    auto const drho_dpCap =
                        (density_plus[d] - density_minus[d]) / (2 * eps_pCap);

                    ASSERT_NEAR(drho_dpCap,
                                density_derivative(d, capillary_pressure),
                                1.e-10);
                }

                // Pertubation of temperature
                auto const eps_T = 0.01;
                auto const T_plus = T + eps_T;

                vars[static_cast<int>(MPL::Variable::capillary_pressure)] =
                    pCap;
                vars[static_cast<int>(MPL::Variable::temperature)] = T_plus;
                vars[static_cast<int>(MPL::Variable::molar_fraction)] =
                    molarFractionFantasy(pGR, pGR_min, pGR_max, pCap, pCap_min,
                                         pCap_max, T_plus, T_min, T_max);

                density_plus = std::get<Eigen::Matrix<double, 3, 1>>(
                    density_model.value(vars, pos, t, dt));

                auto const T_minus = T - eps_T;
                vars[static_cast<int>(MPL::Variable::temperature)] = T_minus;
                vars[static_cast<int>(MPL::Variable::molar_fraction)] =
                    molarFractionFantasy(pGR, pGR_min, pGR_max, pCap, pCap_min,
                                         pCap_max, T_minus, T_min, T_max);

                density_minus = std::get<Eigen::Matrix<double, 3, 1>>(
                    density_model.value(vars, pos, t, dt));

                for (std::size_t d = gas_phase_density; d <= air_density; d++)
                {
                    auto const drho_dT =
                        (density_plus[d] - density_minus[d]) / (2 * eps_T);

                    ASSERT_NEAR(drho_dT, density_derivative(d, temperature),
                                1.e-10);
                }

            }  // end of T-loop
        }      // end of pCap-loop
    }          // end of pGR-loop
}
