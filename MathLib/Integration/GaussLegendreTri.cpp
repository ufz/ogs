// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "GaussLegendreTri.h"

namespace MathLib
{
template <>
const std::array<std::array<double, 2>, GaussLegendreTri<1>::NPoints>
    GaussLegendreTri<1>::X = {{{{1. / 3., 1. / 3.}}}};
template <>
double const GaussLegendreTri<1>::W[1] = {1.0};

const std::array<std::array<double, 2>, GaussLegendreTri<2>::NPoints>
    GaussLegendreTri<2>::X = {
        {{{1. / 6., 1. / 6.}}, {{2. / 3., 1. / 6.}}, {{1. / 6., 2. / 3.}}}};
double const GaussLegendreTri<2>::W[3] = {1. / 3., 1. / 3., 1. / 3.};

const std::array<std::array<double, 2>, GaussLegendreTri<3>::NPoints>
    GaussLegendreTri<3>::X = {{{{1. / 3., 1. / 3.}},
                               {{1. / 5., 3. / 5.}},
                               {{1. / 5., 1. / 5.}},
                               {{3. / 5., 1. / 5.}}}};
double const GaussLegendreTri<3>::W[4] = {-27. / 48., 25. / 48., 25. / 48.,
                                          25. / 48.};

const std::array<std::array<double, 2>, GaussLegendreTri<4>::NPoints>
    GaussLegendreTri<4>::X = {{{{1. / 3., 1. / 3.}},
                               {{0.059715871789770, 0.470142064105115}},
                               {{0.470142064105115, 0.059715871789770}},
                               {{0.470142064105115, 0.470142064105115}},
                               {{0.797426985353087, 0.101286507323456}},
                               {{0.101286507323456, 0.797426985353087}},
                               {{0.101286507323456, 0.101286507323456}}}};
double const GaussLegendreTri<4>::W[7] = {0.225,
                                          0.132394152788506,
                                          0.132394152788506,
                                          0.132394152788506,
                                          0.125939180544827,
                                          0.125939180544827,
                                          0.125939180544827};

}  // namespace MathLib
