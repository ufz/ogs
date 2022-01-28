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

static const double pGR_min = 80000.;
static const double pGR_max = 320000.;
static const double pCap_min = 100.;
static const double pCap_max = 1.e8;
static const double T_min = 270.;
static const double T_max = 520.;

static const double MC = 0.028949;
static const double MW = 0.018052;

double molarMass(const double pGR, const double pCap, const double T)
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

    auto const xnWG =
        a * (pGR - pGR_min) * (pGR - pGR_min) + b * (pGR - pGR_min) +
        c * (pCap - pCap_min) * (pCap - pCap_min) + d * (pCap - pCap_min) +
        e * (T - T_min) * (T - T_min) + f * (T - T_min);
    return xnWG * MW + (1. - xnWG) * MC;
}

std::array<double, 3> molarMassDerivatives(const double pGR,
                                           const double pCap,
                                           const double T)
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

    const double dxn_dpGR = 2 * a * (pGR - pGR_min) + b;
    const double dxn_dpCap = 2 * c * (pCap - pCap_min) + d;
    const double dxn_dT = 2 * e * (T - T_min) + f;

    return {dxn_dpGR * (MW - MC), dxn_dpCap * (MW - MC), dxn_dT * (MW - MC)};
}

TEST(MaterialPropertyLib, IdealGasLawBinaryMixture)
{
    auto const density_model = MPL::IdealGasLawBinaryMixture("density");

    // Enum for dereferencing the derivatives in the return array of the molar
    // mass fantasy function
    enum
    {
        gas_pressure,
        capillary_pressure,
        temperature
    };

    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;

    // Density derivatives of gases are dependent on pressure and temperature.
    // In the case of binary gases, there is also the dependence of the
    // composition (e.g. molar fractions).

    // The average molar mass of the gas phase is represented here by a fantasy
    // function. In Reality, this is actually a physical quantity that depends
    // on the temperature-dependent vapour pressure, the capillary
    // pressure-dependent correction factor for curved menisci and the pressure
    // of the gas phase.

    // In this test, however, the function is merely a polynomial,
    // which over the entire parameter space pGR-pCap-T assumes
    // values in the range [0,1].

    // travel through the entire parameter space
    for (double pGR = pGR_min; pGR <= pGR_max; pGR += 250.)
    {
        for (double pCap = pCap_min; pCap <= pCap_max; pCap *= 2.5)
        {
            for (double T = T_min; T <= T_max; T += 0.25)
            {
                vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR;
                vars[static_cast<int>(MPL::Variable::capillary_pressure)] =
                    pCap;
                vars[static_cast<int>(MPL::Variable::temperature)] = T;

                // Molar mass fantasy function
                auto const MG = molarMass(pGR, pCap, T);

                // ... and its derivatives
                auto const dMG = molarMassDerivatives(pGR, pCap, T);

                vars[static_cast<int>(MPL::Variable::molar_mass)] = MG;

                auto const rhoGR =
                    std::get<double>(density_model.value(vars, pos, t, dt));

                const double R =
                    MaterialLib::PhysicalConstant::IdealGasConstant;

                auto const rhoGR_ = pGR * MG / R / T;
                ASSERT_NEAR(rhoGR, rhoGR_, 1.e-18);

                vars[static_cast<int>(MPL::Variable::molar_mass_derivative)] =
                    dMG[gas_pressure];
                auto const drhoGR_dpGR = std::get<double>(density_model.dValue(
                    vars, MPL::Variable::phase_pressure, pos, t, dt));

                vars[static_cast<int>(MPL::Variable::molar_mass_derivative)] =
                    dMG[capillary_pressure];
                auto const drhoGR_dpCap = std::get<double>(density_model.dValue(
                    vars, MPL::Variable::capillary_pressure, pos, t, dt));

                vars[static_cast<int>(MPL::Variable::molar_mass_derivative)] =
                    dMG[temperature];
                auto const drhoGR_dT = std::get<double>(density_model.dValue(
                    vars, MPL::Variable::temperature, pos, t, dt));

                // Check derivatives via central differences
                auto central_difference =
                    [&](double plus, double minus, double eps)
                    {
                        return (plus - minus) / (2 * eps);
                    };

                // Pertubation of gas pressure
                auto const eps_pGR = 10.;

                vars[static_cast<int>(MPL::Variable::phase_pressure)] =
                    pGR + eps_pGR;
                vars[static_cast<int>(MPL::Variable::molar_mass)] =
                    molarMass(pGR + eps_pGR, pCap, T);

                auto rhoGR_plus =
                    std::get<double>(density_model.value(vars, pos, t, dt));

                vars[static_cast<int>(MPL::Variable::phase_pressure)] =
                    pGR - eps_pGR;
                vars[static_cast<int>(MPL::Variable::molar_mass)] =
                    molarMass(pGR - eps_pGR, pCap, T);

                auto rhoGR_minus =
                    std::get<double>(density_model.value(vars, pos, t, dt));

                ASSERT_NEAR(
                    drhoGR_dpGR,
                    central_difference(rhoGR_plus, rhoGR_minus, eps_pGR),
                    1.e-10);

                vars[static_cast<int>(MPL::Variable::phase_pressure)] = pGR;

                // Pertubation of capillary pressure
                auto const eps_pCap = 1.;

                vars[static_cast<int>(MPL::Variable::capillary_pressure)] =
                    pCap + eps_pCap;
                vars[static_cast<int>(MPL::Variable::molar_mass)] =
                    molarMass(pGR, pCap + eps_pCap, T);

                rhoGR_plus =
                    std::get<double>(density_model.value(vars, pos, t, dt));

                vars[static_cast<int>(MPL::Variable::capillary_pressure)] =
                    pCap - eps_pCap;
                vars[static_cast<int>(MPL::Variable::molar_mass)] =
                    molarMass(pGR, pCap - eps_pCap, T);

                rhoGR_minus =
                    std::get<double>(density_model.value(vars, pos, t, dt));

                ASSERT_NEAR(
                    drhoGR_dpCap,
                    central_difference(rhoGR_plus, rhoGR_minus, eps_pCap),
                    1.e-10);

                vars[static_cast<int>(MPL::Variable::capillary_pressure)] =
                    pCap;

                // Pertubation of temperature
                auto const eps_T = .01;

                vars[static_cast<int>(MPL::Variable::temperature)] = T + eps_T;
                vars[static_cast<int>(MPL::Variable::molar_mass)] =
                    molarMass(pGR, pCap, T + eps_T);

                rhoGR_plus =
                    std::get<double>(density_model.value(vars, pos, t, dt));

                vars[static_cast<int>(MPL::Variable::temperature)] = T - eps_T;
                vars[static_cast<int>(MPL::Variable::molar_mass)] =
                    molarMass(pGR, pCap, T - eps_T);

                rhoGR_minus =
                    std::get<double>(density_model.value(vars, pos, t, dt));

                ASSERT_NEAR(drhoGR_dT,
                            central_difference(rhoGR_plus, rhoGR_minus, eps_T),
                            1.e-10);

            }  // end of T-loop
        }      // end of pCap-loop
    }          // end of pGR-loop
}
