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

#include "ExtrapolationGaussQuad.h"

#include <cmath>
#include <limits>
#include <stdexcept>

#include "MathLib/Integration/GaussLegendre.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/Integration/IntegrationGaussRegular.h"

namespace NumLib
{

std::size_t ExtrapolationGaussQuad::getNodeIndexOfGaussQuad(std::size_t nGaussLevel, std::size_t igp)
{
    std::size_t LoIndex = -1;
    static const double epsilon = std::numeric_limits<double>::epsilon();
    double r[3] = {.0};
    IntegrationGaussRegular<2>::getPoint(nGaussLevel, igp, r);

    if (r[0] > 0.0 && r[1] > 0.0)
        LoIndex = 0;
    else if (r[0] < 0.0 && r[1] > 0.0)
        LoIndex = 1;
    else if (r[0] < 0.0 && r[1] < 0.0)
        LoIndex = 2;
    else if (r[0] > 0.0 && r[1] < 0.0)
        LoIndex = 3;
    else if (std::abs(r[0]) < epsilon && r[1] > 0.0)
        LoIndex = 4;
    else if (r[0] < 0.0 && std::abs(r[1]) < epsilon)
        LoIndex = 5;
    else if (std::abs(r[0]) < epsilon && r[1] < 0.0)
        LoIndex = 6;
    else if (r[0] > 0.0 && std::abs(r[1]) < epsilon)
        LoIndex = 7;
    else if (std::abs(r[0]) < epsilon && std::abs(r[1]) < epsilon)
        LoIndex = 8;

    return LoIndex;
}

double ExtrapolationGaussQuad::calculateXi_p(std::size_t nGaussLevel)
{
    double Xi_p = 0.0;
    for (std::size_t gp=0; gp<nGaussLevel; gp++)
        Xi_p = std::max(Xi_p, std::abs(MathLib::GaussLegendre::getPoint(nGaussLevel, gp).first));
    return 1.0 / Xi_p; // 1 means nodes
}

void ExtrapolationGaussQuad::getExtrapolatedPoints(std::size_t nodeIdOfGaussQuad, double Xi_p, double* r)
{
    switch (nodeIdOfGaussQuad)
    {
        case 0:
            r[0] = Xi_p;
            r[1] = Xi_p;
            break;
        case 1:
            r[0] = -Xi_p;
            r[1] = Xi_p;
            break;
        case 2:
            r[0] = -Xi_p;
            r[1] = -Xi_p;
            break;
        case 3:
            r[0] = Xi_p;
            r[1] = -Xi_p;
            break;
        default:
            throw std::logic_error("ExtrapolationGaussQuad: only four nodes are supported.");
    }
}


} //end namespace

