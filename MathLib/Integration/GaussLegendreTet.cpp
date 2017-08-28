/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GaussLegendreTet.h"

namespace MathLib
{
template <>
const std::array<std::array<double, 3>, GaussLegendreTet<1>::NPoints>
    GaussLegendreTet<1>::X = {{{{1. / 4., 1. / 4., 1. / 4.}}}};
template <>
double const GaussLegendreTet<1>::W[1] = {1. / 6.};

const std::array<std::array<double, 3>, GaussLegendreTet<2>::NPoints>
GaussLegendreTet<2>::X = {{ {{1./4., 1./4., 1./4.}},
                            {{1./6., 1./6., 1./6.}},
                            {{1./2., 1./6., 1./6.}},
                            {{1./6., 1./2., 1./6.}},
                            {{1./6., 1./6., 1./2.}} }};
double const GaussLegendreTet<2>::W[5] = {-2./15., 0.075, 0.075, 0.075, 0.075};

std::array<std::array<double, 3>, GaussLegendreTet<3>::NPoints> initGLTet3X()
{
    // Cf. Gellert, M., Harbord, R., 1991. Moderate degree cubature formulas for
    // 3-D tetrahedral finite-element approximations. Communications in Applied
    // Numerical Methods 7, 487–495. doi:10.1002/cnm.1630070609
    const double a = 0.0673422422100983;
    const double b = 0.3108859192633005;
    const double c = 0.7217942490673264;
    const double d = 0.0927352503108912;
    const double e = 0.4544962958743506;
    const double f = 0.045503704125649;

    return {{{a, b, b},
             {b, a, b},
             {b, b, a},
             {b, b, b},
             {c, d, d},
             {d, c, d},
             {d, d, c},
             {d, d, d},
             {e, e, f},
             {e, f, e},
             {e, f, f},
             {f, e, e},
             {f, e, f},
             {f, f, e}}};
}

const std::array<std::array<double, 3>, GaussLegendreTet<3>::NPoints>
    GaussLegendreTet<3>::X = initGLTet3X();

static const double p = 0.1126879257180162 / 6.;
static const double q = 0.0734930431163619 / 6.;
static const double r = 0.0425460207770812 / 6.;

double const GaussLegendreTet<3>::W[GaussLegendreTet<3>::NPoints] = {
    p, p, p, p, q, q, q, q, r, r, r, r, r, r};
}
