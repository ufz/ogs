/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 4, 2021, 3:05 PM
 */

#include "WaterThermalConductivityIAPWS.h"

#include <cmath>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
// Li and Lij are coefficients from Tables 1 and 2 (Daucik and Dooley, 2011)
static const double Li[5] = {2.443221e-3, 1.323095e-2, 6.770357e-3,
                             -3.454586e-3, 4.096266e-4};
static const double Lij[5][6] = {
    {1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634,
     0.00609859258},
    {2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.0832827019,
     -0.00719201245},
    {2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278,
     -0.0205938816},
    {-1.21051378, 1.60812989, -0.621178141, 0.0716373224, 0, 0},
    {-2.7203370, 4.57586331, -3.18369245, 1.1168348, -0.19268305, 0.012913842}};

static double computeBarLambda0Factor(const double barT);

static std::array<double, 5> computeSeriesFactorTForLambda1(const double barT);
static std::array<double, 6> computeSeriesFactorRhoForLambda1(
    const double bar_rho);
static double computeBarLambda1Factor(
    const std::array<double, 5>& series_factorT,
    const std::array<double, 6>& series_factorRho);

static double computedBarLambda_dbarT(const double barT, double bar_rho);
static double computedBarLambda_dbarRho(const double barT, double bar_rho);

PropertyDataType WaterThermalConductivityIAPWS::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double bar_T = variable_array.temperature / ref_T_;
    const double bar_rho = variable_array.density / ref_rho_;

    const double lambda0 = std::sqrt(bar_T) / computeBarLambda0Factor(bar_T);

    const auto& series_factorT = computeSeriesFactorTForLambda1(bar_T);
    const auto& series_factorRho = computeSeriesFactorRhoForLambda1(bar_rho);
    const double lambda1 = std::exp(
        bar_rho * computeBarLambda1Factor(series_factorT, series_factorRho));

    return lambda0 * lambda1 * ref_lambda_;
}

PropertyDataType WaterThermalConductivityIAPWS::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double bar_T = variable_array.temperature / ref_T_;
    const double bar_rho = variable_array.density / ref_rho_;

    switch (variable)
    {
        case Variable::temperature:
            return ref_lambda_ * computedBarLambda_dbarT(bar_T, bar_rho) /
                   ref_T_;
        case Variable::density:
            return ref_lambda_ * computedBarLambda_dbarRho(bar_T, bar_rho) /
                   ref_rho_;
        default:
            OGS_FATAL(
                "WaterThermalConductivityIAPWS::dValue is implemented for "
                "derivatives with respect to temperature and liquid density "
                "only.");
    }
}

double computeBarLambda0Factor(const double barT)
{
    double sum_val = 0.;
    double barT_i = 1.;
    for (double value : Li)
    {
        sum_val += (value / barT_i);
        barT_i *= barT;
    }
    return sum_val;
}

std::array<double, 5> computeSeriesFactorTForLambda1(const double barT)
{
    std::array<double, 5> series_factorT;
    series_factorT[0] = 1.;
    const double barT_fac = 1 / barT - 1.0;
    for (int i = 1; i < 5; i++)
    {
        series_factorT[i] = series_factorT[i - 1] * barT_fac;
    }

    return series_factorT;
}

std::array<double, 6> computeSeriesFactorRhoForLambda1(const double bar_rho)
{
    std::array<double, 6> series_factorRho;
    series_factorRho[0] = 1.;
    for (int i = 1; i < 6; i++)
    {
        series_factorRho[i] = series_factorRho[i - 1] * (bar_rho - 1.0);
    }

    return series_factorRho;
}

double computeBarLambda1Factor(const std::array<double, 5>& series_factorT,
                               const std::array<double, 6>& series_factorRho)
{
    double sum_val = 0.;
    for (int i = 0; i < 5; i++)
    {
        double sum_val_j = 0;
        for (int j = 0; j < 6; j++)
        {
            sum_val_j += Lij[i][j] * series_factorRho[j];
        }
        sum_val += series_factorT[i] * sum_val_j;
    }

    return sum_val;
}

double computedBarLambda_dbarT(const double barT, double bar_rho)
{
    const double lambda0_factor = computeBarLambda0Factor(barT);
    const double sqrt_barT = std::sqrt(barT);

    double dlambda0_factor_dbarT = 0.0;
    double barT_i = barT * barT;
    for (int i = 1; i < 5; i++)
    {
        dlambda0_factor_dbarT -= static_cast<double>(i) * (Li[i] / barT_i);
        barT_i *= barT;
    }

    const double dbar_lambda0_dbarT =
        0.5 / (lambda0_factor * sqrt_barT) -
        sqrt_barT * dlambda0_factor_dbarT / (lambda0_factor * lambda0_factor);

    const auto& series_factorT = computeSeriesFactorTForLambda1(barT);
    const auto& series_factorRho = computeSeriesFactorRhoForLambda1(bar_rho);

    double dlambda1_factor_dbarT = 0.0;
    for (int i = 1; i < 5; i++)
    {
        double sum_val_j = 0;
        for (int j = 0; j < 6; j++)
        {
            sum_val_j += Lij[i][j] * series_factorRho[j];
        }
        dlambda1_factor_dbarT -= static_cast<double>(i) *
                                 series_factorT[i - 1] * sum_val_j /
                                 (barT * barT);
    }

    const double lambda1_factor =
        computeBarLambda1Factor(series_factorT, series_factorRho);
    const double dbar_lambda1_dbarT =
        bar_rho * std::exp(bar_rho * lambda1_factor) * dlambda1_factor_dbarT;

    return dbar_lambda0_dbarT * std::exp(bar_rho * lambda1_factor) +
           dbar_lambda1_dbarT * sqrt_barT / lambda0_factor;
}

double computedBarLambda_dbarRho(const double barT, double bar_rho)
{
    const auto& series_factorT = computeSeriesFactorTForLambda1(barT);
    const auto& series_factorRho = computeSeriesFactorRhoForLambda1(bar_rho);

    double dlambda1_factor_dbar_rho = 0.0;
    for (int i = 0; i < 5; i++)
    {
        double sum_val_j = 0;
        for (int j = 1; j < 6; j++)
        {
            sum_val_j +=
                static_cast<double>(j) * Lij[i][j] * series_factorRho[j - 1];
        }
        dlambda1_factor_dbar_rho += series_factorT[i] * sum_val_j;
    }

    const double lambda0 = std::sqrt(barT) / computeBarLambda0Factor(barT);

    const double lambda1_factor =
        computeBarLambda1Factor(series_factorT, series_factorRho);
    return lambda0 * std::exp(bar_rho * lambda1_factor) *
           (lambda1_factor + bar_rho * dlambda1_factor_dbar_rho);
}

}  // namespace MaterialPropertyLib
