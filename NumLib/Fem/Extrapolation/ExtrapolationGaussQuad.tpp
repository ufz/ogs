/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-08-13
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include <cassert>
#include <cmath>
#include <stdexcept>
#include "MathLib/Integration/GaussLegendre.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/Integration/IntegrationGaussRegular.h"

namespace NumLib
{

template <class T_GP_VEC, class T_NOD_VEC>
void ExtrapolationGaussQuad::extrapolate(const T_GP_VEC &gp_values, T_NOD_VEC &nodal_values)
{
    static const std::size_t nExtrapolatedNodes = 4; // only corner nodes
    const std::size_t nGaussPoints = gp_values.size();
    assert(nGaussPoints >= nExtrapolatedNodes);
    const std::size_t nGaussLevel = std::sqrt(nGaussPoints);

    // reorder gauss point values
    T_NOD_VEC reordered_gp_values(nExtrapolatedNodes);
    for (int i=0; i<gp_values.size(); i++) {
        std::size_t nod_id = getNodeIndexOfGaussQuad(nGaussLevel, i);
        if (nod_id < nExtrapolatedNodes)
            reordered_gp_values[nod_id] = gp_values[i];
    }

    // calculate Xi_p
    const double Xi_p = calculateXi_p(nGaussLevel);

    // extrapolate linearly
    nodal_values.resize(nExtrapolatedNodes);
    double x[3] = {};
    T_NOD_VEC N(nExtrapolatedNodes);
    for (std::size_t i=0; i<nExtrapolatedNodes; i++) {
        getExtrapolatedPoints(i, Xi_p, x);
        ShapeQuad4::computeShapeFunction(x, N);
        nodal_values[i] = N.dot(reordered_gp_values);
    }
}

} // end namespace

