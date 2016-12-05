/**
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   WaterViscosityIAPWS.cpp
 *
 * Created on December 1, 2016, 1:41 PM
 */

#include "WaterViscosityIAPWS.h"

#include <cmath>

namespace MaterialLib
{
namespace Fluid
{
double WaterViscosityIAPWS::getValue(const ArrayType& var_vals) const
{
    const double bar_T =
        var_vals[static_cast<unsigned>(PropertyVariableType::T)] / _ref_T;
    const double bar_rho =
        var_vals[static_cast<unsigned>(PropertyVariableType::rho)] / _ref_rho;

    const double mu0 = 100. * std::sqrt(bar_T) / computeBarMu0Factor(bar_T);

    computeSeriesFactorForMu1(bar_T, bar_rho);
    const double mu1 = std::exp(bar_rho * computeBarMu1Factor());

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
            computeSeriesFactorForMu1(bar_T, bar_rho);
            return _ref_mu * computedBarMu_dbarRho(bar_T, bar_rho) / _ref_rho;
        default:
            return 0.;
    }
}

double WaterViscosityIAPWS::computeBarMu0Factor(const double barT) const
{
    double sum_val = 0.;
    double barT_i = 1.;
    for (int i = 0; i < 4; i++)
    {
        sum_val += (_hi[i] / barT_i);
        barT_i *= barT;
    }
    return sum_val;
}

void WaterViscosityIAPWS::computeSeriesFactorForMu1(const double barT,
                                                    const double bar_rho) const
{
    _series_factorT[0] = 1.;
    const double barT_fac = 1 / barT - 1.0;
    for (int i = 1; i < 6; i++)
    {
        _series_factorT[i] = _series_factorT[i - 1] * barT_fac;
    }

    _series_factorRho[0] = 1.;
    for (int i = 1; i < 7; i++)
    {
        _series_factorRho[i] = _series_factorRho[i - 1] * (bar_rho - 1.0);
    }
}

double WaterViscosityIAPWS::computeBarMu1Factor() const
{
    double sum_val = 0.;
    for (int i = 0; i < 6; i++)
    {
        double sum_val_j = 0;
        for (int j = 0; j < 7; j++)
        {
            sum_val_j += _hij[i][j] * _series_factorRho[j];
        }
        sum_val += _series_factorT[i] * sum_val_j;
    }

    return sum_val;
}

double WaterViscosityIAPWS::computedBarMu_dbarT(const double barT,
                                                double bar_rho) const
{
    const double mu0_factor = computeBarMu0Factor(barT);
    const double sqrt_barT = std::sqrt(barT);

    double dmu0_factor_dbarT = 0.0;
    double barT_i = barT * barT;
    for (int i = 1; i < 4; i++)
    {
        dmu0_factor_dbarT -= static_cast<double>(i) * (_hi[i] / barT_i);
        barT_i *= barT;
    }

    const double dbar_mu0_dbarT =
        50. / (mu0_factor * sqrt_barT) -
        100. * sqrt_barT * dmu0_factor_dbarT / (mu0_factor * mu0_factor);

    double dmu1_factor_dbarT = 0.0;
    for (int i = 1; i < 6; i++)
    {
        double sum_val_j = 0;
        for (int j = 0; j < 7; j++)
        {
            sum_val_j += _hij[i][j] * _series_factorRho[j];
        }
        dmu1_factor_dbarT -= static_cast<double>(i) * _series_factorT[i - 1] *
                             sum_val_j / (barT * barT);
    }

    const double mu1_factor = computeBarMu1Factor();
    const double dbar_mu1_dbarT =
        bar_rho * std::exp(bar_rho * mu1_factor) * dmu1_factor_dbarT;

    return dbar_mu0_dbarT * std::exp(bar_rho * mu1_factor) +
           dbar_mu1_dbarT * 100. * sqrt_barT / mu0_factor;
}

double WaterViscosityIAPWS::computedBarMu_dbarRho(const double barT,
                                                  double bar_rho) const
{
    double dmu1_factor_dbar_rho = 0.0;
    for (int i = 0; i < 6; i++)
    {
        double sum_val_j = 0;
        for (int j = 1; j < 7; j++)
        {
            sum_val_j +=
                static_cast<double>(j) * _hij[i][j] * _series_factorRho[j - 1];
        }
        dmu1_factor_dbar_rho += _series_factorT[i] * sum_val_j;
    }

    const double mu0 = 100. * std::sqrt(barT) / computeBarMu0Factor(barT);

    const double mu1_factor = computeBarMu1Factor();
    return mu0 * std::exp(bar_rho * mu1_factor) *
           (mu1_factor + bar_rho * dmu1_factor_dbar_rho);
}

}  // end namespace
}  // end namespace
