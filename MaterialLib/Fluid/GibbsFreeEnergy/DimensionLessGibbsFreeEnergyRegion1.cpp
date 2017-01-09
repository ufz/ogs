/**
 *  \brief Define members of the class for dimensionless Gibbs free energy.
 *
 *  \copyright
 *   Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   DimensionLessGibbsFreeEnergyRegion1.cpp
 *
 * Created on December 8, 2016, 12:31 PM
 */

#include "DimensionLessGibbsFreeEnergyRegion1.h"

#include <cmath>

namespace MaterialLib
{
namespace Fluid
{
static const double ni[34] = {
    1.4632971213167E-01,  -8.4548187169114E-01, -3.7563603672040E+00,
    3.3855169168385E+00,  -9.5791963387872E-01, 1.5772038513228E-01,
    -1.6616417199501E-02, 8.1214629983568E-04,  2.8319080123804E-04,
    -6.0706301565874E-04, -1.8990068218419E-02, -3.2529748770505E-02,
    -2.1841717175414E-02, -5.2838357969930E-05, -4.7184321073267E-04,
    -3.0001780793026E-04, 4.7661393906987E-05,  -4.4141845330846E-06,
    -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06,
    -8.5205128120103E-10, -2.2425281908000E-06, -6.5171222895601E-07,
    -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09,
    -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20,
    2.6335781662795E-23,  -1.1947622640071E-23, 1.8228094581404E-24,
    -9.35370872924580E-26};

static const double li[34] = {0, 0, 0, 0, 0,  0,  0,  0,  1,  1, 1, 1,
                              1, 1, 2, 2, 2,  2,  2,  3,  3,  3, 4, 4,
                              4, 5, 8, 8, 21, 23, 29, 30, 31, 32};

static const double ji[34] = {
    -2, -1, 0,  1, 2, 3,  4,  5,  -9, -7,  -1, 0,   1,   3,   -3,  0,   1,
    3,  17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41};

double DimensionLessGibbsFreeEnergyRegion1::get_gamma(const double tau,
                                                      const double pi) const
{
    double val = 0.;
    for (int i = 0; i < 34; i++)
    {
        val += ni[i] * std::pow(7.1 - pi, li[i]) * std::pow(tau - 1.222, ji[i]);
    }

    return val;
}

double DimensionLessGibbsFreeEnergyRegion1::get_dgamma_dtau(
    const double tau, const double pi) const
{
    double val = 0.;
    for (int i = 0; i < 34; i++)
    {
        val += ni[i] * ji[i] * std::pow(7.1 - pi, li[i]) *
               std::pow(tau - 1.222, ji[i] - 1.0);
    }

    return val;
}

double DimensionLessGibbsFreeEnergyRegion1::get_dgamma_dtau_dtau(
    const double tau, const double pi) const
{
    double val = 0.;
    for (int i = 0; i < 34; i++)
    {
        val += ni[i] * ji[i] * (ji[i] - 1.0) * std::pow(7.1 - pi, li[i]) *
               std::pow(tau - 1.222, ji[i] - 2.0);
    }

    return val;
}

double DimensionLessGibbsFreeEnergyRegion1::get_dgamma_dpi(
    const double tau, const double pi) const
{
    double val = 0.;
    for (int i = 0; i < 34; i++)
    {
        val -= ni[i] * li[i] * std::pow(7.1 - pi, li[i] - 1.0) *
               std::pow(tau - 1.222, ji[i]);
    }

    return val;
}

double DimensionLessGibbsFreeEnergyRegion1::get_dgamma_dpi_dpi(
    const double tau, const double pi) const
{
    double val = 0.;
    for (int i = 0; i < 34; i++)
    {
        val += ni[i] * li[i] * (li[i] - 1.0) * std::pow(7.1 - pi, li[i] - 2.0) *
               std::pow(tau - 1.222, ji[i]);
    }

    return val;
}

double DimensionLessGibbsFreeEnergyRegion1::get_dgamma_dtau_dpi(
    const double tau, const double pi) const
{
    double val = 0.;
    for (int i = 0; i < 34; i++)
    {
        val -= ni[i] * ji[i] * li[i] * std::pow(7.1 - pi, li[i] - 1.0) *
               std::pow(tau - 1.222, ji[i] - 1.0);
    }

    return val;
}

}  // end namespace
}  // end namespace
