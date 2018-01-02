/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-08-13
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GaussLegendreTri.h"

namespace MathLib
{

template <> const std::array<std::array<double, 2>, GaussLegendreTri<1>::NPoints> GaussLegendreTri<1>::X = {{{{1./3., 1./3.}}}};
template <> double const GaussLegendreTri<1>::W[1] = {1.0};

const std::array<std::array<double, 2>, GaussLegendreTri<2>::NPoints> GaussLegendreTri<2>::X = {{{{1./6., 1./6.}}, {{2./3., 1./6.}}, {{1./6., 2./3.}}}};
double const GaussLegendreTri<2>::W[3] = {1./3., 1./3., 1./3.};

const std::array<std::array<double, 2>, GaussLegendreTri<3>::NPoints> GaussLegendreTri<3>::X = {{{{1./3., 1./3.}}, {{1./5., 3./5.}}, {{1./5., 1./5.}}, {{3./5., 1./5.}}}};
double const GaussLegendreTri<3>::W[4] = {-27./48., 25./48., 25./48., 25./48.};

}
