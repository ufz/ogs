/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
    GaussLegendreTet<2>::X = {{{{1. / 4., 1. / 4., 1. / 4.}},
                               {{1. / 6., 1. / 6., 1. / 6.}},
                               {{1. / 2., 1. / 6., 1. / 6.}},
                               {{1. / 6., 1. / 2., 1. / 6.}},
                               {{1. / 6., 1. / 6., 1. / 2.}}}};
double const GaussLegendreTet<2>::W[5] = {-2. / 15., 0.075, 0.075, 0.075,
                                          0.075};

static std::array<std::array<double, 3>, GaussLegendreTet<3>::NPoints>
initGLTet3X()
{
    // Cf. Gellert, M., Harbord, R., 1991. Moderate degree cubature formulas for
    // 3-D tetrahedral finite-element approximations. Communications in Applied
    // Numerical Methods 7, 487-495. doi:10.1002/cnm.1630070609
    const double a = 0.0673422422100983;
    const double b = 0.3108859192633005;
    const double c = 0.7217942490673264;
    const double d = 0.0927352503108912;
    const double e = 0.4544962958743506;
    const double f = 0.045503704125649;

    return {{{{a, b, b}},
             {{b, a, b}},
             {{b, b, a}},
             {{b, b, b}},
             {{c, d, d}},
             {{d, c, d}},
             {{d, d, c}},
             {{d, d, d}},
             {{e, e, f}},
             {{e, f, e}},
             {{e, f, f}},
             {{f, e, e}},
             {{f, e, f}},
             {{f, f, e}}}};
}

const std::array<std::array<double, 3>, GaussLegendreTet<3>::NPoints>
    GaussLegendreTet<3>::X = initGLTet3X();

static const double p = 0.1126879257180162 / 6.;
static const double q = 0.0734930431163619 / 6.;
static const double r = 0.0425460207770812 / 6.;

double const GaussLegendreTet<3>::W[GaussLegendreTet<3>::NPoints] = {
    p, p, p, p, q, q, q, q, r, r, r, r, r, r};

static std::array<std::array<double, 3>, GaussLegendreTet<4>::NPoints>
initGLTet4X()
{
    // Cf. Gellert, M., Harbord, R., 1991. Moderate degree cubature formulas for
    // 3-D tetrahedral finite-element approximations. Communications in Applied
    // Numerical Methods 7, 487-495. doi:10.1002/cnm.1630070609
    const double a = 0.3797582452067875;
    const double b = 0.1202417547932126;

    return {{{{0.0, 1. / 3, 1. / 3}},
             {{1. / 3, 0.0, 1. / 3}},
             {{1. / 3, 1. / 3, 0.0}},
             {{1. / 3, 1. / 3, 1. / 3}},
             {{8. / 11, 1. / 11, 1. / 11}},
             {{1. / 11, 8. / 11, 1. / 11}},
             {{1. / 11, 1. / 11, 8. / 11}},
             {{1. / 11, 1. / 11, 1. / 11}},
             {{0.0, 0.0, 0.5}},
             {{0.0, 0.5, 0.0}},
             {{0.0, 0.5, 0.5}},
             {{0.5, 0.0, 0.0}},
             {{0.5, 0.0, 0.5}},
             {{0.5, 0.5, 0.0}},
             {{a, a, b}},
             {{a, b, a}},
             {{a, b, b}},
             {{b, a, a}},
             {{b, a, b}},
             {{b, b, a}}}};
}
const std::array<std::array<double, 3>, GaussLegendreTet<4>::NPoints>
    GaussLegendreTet<4>::X = initGLTet4X();

static const double s = 81. / 2240. / 6.;
static const double t = 161051. / 2304960. / 6.;
static const double u = 409. / 31395. / 6.;
static const double v = 2679769. / 32305455. / 6.;

double const GaussLegendreTet<4>::W[GaussLegendreTet<4>::NPoints] = {
    s, s, s, s, t, t, t, t, u, u, u, u, u, u, v, v, v, v, v, v};

}  // namespace MathLib
