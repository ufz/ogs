/**
 *  \copyright
 *   Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   WaterViscosityIAPWS.cpp
 *
 * Created on December 1, 2016, 1:41 PM
 */

#include "WaterViscosityIAPWS.h"

#include <array>
#include <cmath>

namespace MaterialLib
{
namespace Fluid
{
static const double Hi[4] = {1.67752, 2.20462, 0.6366564, -0.241605};
static const double Hij[6][7] = {
    {0.520094, 0.222531, -0.281378, 0.161913, -0.0325372, 0, 0},
    {0.0850895, 0.999115, -0.906851, 0.257399, 0, 0, 0},
    {-1.08374, 1.88797, -0.772479, 0, 0, 0, 0},
    {-0.289555, 1.26613, -0.489837, 0, 0.0698452, 0, -0.00435673},
    {0, 0, -0.25704, 0, 0, 0.00872102, 0},
    {0, 0.120573, 0, 0, 0, 0, -0.000593264}};

static double computeBarMu0Factor(const double barT);

static std::array<double, 6> computeSeriesFactorTForMu1(const double barT);
static std::array<double, 7> computeSeriesFactorRhoForMu1(const double bar_rho);
static double computeBarMu1Factor(
    const std::array<double, 6>& series_factorT,
    const std::array<double, 7>& series_factorRho);

static double computedBarMu_dbarT(const double barT, double bar_rho);
static double computedBarMu_dbarRho(const double barT, double bar_rho);

double WaterViscosityIAPWS::getValue(const ArrayType& var_vals) const
{
    const double bar_T =
        var_vals[static_cast<unsigned>(PropertyVariableType::T)] / _ref_T;
    const double bar_rho =
        var_vals[static_cast<unsigned>(PropertyVariableType::rho)] / _ref_rho;

    const double mu0 = 100. * std::sqrt(bar_T) / computeBarMu0Factor(bar_T);

    const auto& series_factorT = computeSeriesFactorTForMu1(bar_T);
    const auto& series_factorRho = computeSeriesFactorRhoForMu1(bar_rho);
    const double mu1 = std::exp(
        bar_rho * computeBarMu1Factor(series_factorT, series_factorRho));

    return mu0 * mu1 * _ref_mu;
}

double WaterViscosityIAPWS::getdValue(const ArrayType& var_vals,
                                      const PropertyVariableType var_type) const
{
    const double bar_T =
        var_vals[static_cast<unsigned>(PropertyVariableType::T)] / _ref_T;
    const double bar_rho =
        var_vals[static_cast<unsigned>(PropertyVariableType::rho)] / _ref_rho;

    switch (var_type)
    {
        case PropertyVariableType::T:
            return _ref_mu * computedBarMu_dbarT(bar_T, bar_rho) / _ref_T;
        case PropertyVariableType::rho:
            return _ref_mu * computedBarMu_dbarRho(bar_T, bar_rho) / _ref_rho;
        default:
            return 0.;
    }
}

double computeBarMu0Factor(const double barT)
{
    double sum_val = 0.;
    double barT_i = 1.;
    for (double i : Hi)
    {
        sum_val += (i / barT_i);
        barT_i *= barT;
    }
    return sum_val;
}

std::array<double, 6> computeSeriesFactorTForMu1(const double barT)
{
    std::array<double, 6> series_factorT;
    series_factorT[0] = 1.;
    const double barT_fac = 1 / barT - 1.0;
    for (int i = 1; i < 6; i++)
    {
        series_factorT[i] = series_factorT[i - 1] * barT_fac;
    }

    return series_factorT;
}

std::array<double, 7> computeSeriesFactorRhoForMu1(const double bar_rho)
{
    std::array<double, 7> series_factorRho;
    series_factorRho[0] = 1.;
    for (int i = 1; i < 7; i++)
    {
        series_factorRho[i] = series_factorRho[i - 1] * (bar_rho - 1.0);
    }

    return series_factorRho;
}

double computeBarMu1Factor(const std::array<double, 6>& series_factorT,
                           const std::array<double, 7>& series_factorRho)
{
    double sum_val = 0.;
    for (int i = 0; i < 6; i++)
    {
        double sum_val_j = 0;
        for (int j = 0; j < 7; j++)
        {
            sum_val_j += Hij[i][j] * series_factorRho[j];
        }
        sum_val += series_factorT[i] * sum_val_j;
    }

    return sum_val;
}

double computedBarMu_dbarT(const double barT, double bar_rho)
{
    const double mu0_factor = computeBarMu0Factor(barT);
    const double sqrt_barT = std::sqrt(barT);

    double dmu0_factor_dbarT = 0.0;
    double barT_i = barT * barT;
    for (int i = 1; i < 4; i++)
    {
        dmu0_factor_dbarT -= static_cast<double>(i) * (Hi[i] / barT_i);
        barT_i *= barT;
    }

    const double dbar_mu0_dbarT =
        50. / (mu0_factor * sqrt_barT) -
        100. * sqrt_barT * dmu0_factor_dbarT / (mu0_factor * mu0_factor);

    const auto& series_factorT = computeSeriesFactorTForMu1(barT);
    const auto& series_factorRho = computeSeriesFactorRhoForMu1(bar_rho);

    double dmu1_factor_dbarT = 0.0;
    for (int i = 1; i < 6; i++)
    {
        double sum_val_j = 0;
        for (int j = 0; j < 7; j++)
        {
            sum_val_j += Hij[i][j] * series_factorRho[j];
        }
        dmu1_factor_dbarT -= static_cast<double>(i) * series_factorT[i - 1] *
                             sum_val_j / (barT * barT);
    }

    const double mu1_factor =
        computeBarMu1Factor(series_factorT, series_factorRho);
    const double dbar_mu1_dbarT =
        bar_rho * std::exp(bar_rho * mu1_factor) * dmu1_factor_dbarT;

    return dbar_mu0_dbarT * std::exp(bar_rho * mu1_factor) +
           dbar_mu1_dbarT * 100. * sqrt_barT / mu0_factor;
}

double computedBarMu_dbarRho(const double barT, double bar_rho)
{
    const auto& series_factorT = computeSeriesFactorTForMu1(barT);
    const auto& series_factorRho = computeSeriesFactorRhoForMu1(bar_rho);

    double dmu1_factor_dbar_rho = 0.0;
    for (int i = 0; i < 6; i++)
    {
        double sum_val_j = 0;
        for (int j = 1; j < 7; j++)
        {
            sum_val_j +=
                static_cast<double>(j) * Hij[i][j] * series_factorRho[j - 1];
        }
        dmu1_factor_dbar_rho += series_factorT[i] * sum_val_j;
    }

    const double mu0 = 100. * std::sqrt(barT) / computeBarMu0Factor(barT);

    const double mu1_factor =
        computeBarMu1Factor(series_factorT, series_factorRho);
    return mu0 * std::exp(bar_rho * mu1_factor) *
           (mu1_factor + bar_rho * dmu1_factor_dbar_rho);
}

}  // end namespace
}  // end namespace
