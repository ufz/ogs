/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include "GaussLegendre.h"

namespace MathLib
{

template <> double const GaussLegendre<1>::X[1] = {0.};
template <> double const GaussLegendre<1>::W[1] = {2.};

template <> double const GaussLegendre<2>::X[2] = {0.577350269189626, -0.577350269189626};
template <> double const GaussLegendre<2>::W[2] = {1., 1.};

template <> double const GaussLegendre<3>::X[3] = {0.774596669241483, 0., -0.774596669241483};
template <> double const GaussLegendre<3>::W[3] = {5./9, 8./9, 5./9};

template <> double const GaussLegendre<4>::X[4] =
        {-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053};
template <> double const GaussLegendre<4>::W[4] =
        { 0.347854845137454,  0.652145154862546, 0.652145154862546, 0.347854845137454};

}   // namespace MathLib
