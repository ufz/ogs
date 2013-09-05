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


#include "MathLib/Integration/GaussLegendre.h"

namespace NumLib
{

template <>
inline std::tuple<std::size_t, std::size_t, std::size_t>
IntegrationGaussRegular<2>::getPosition(std::size_t nGauss, std::size_t igp)
{
    std::size_t gp_r = igp / nGauss;
    std::size_t gp_s = igp % nGauss;
    return std::make_tuple(gp_r, gp_s, 0u);
}

template <>
inline std::tuple<std::size_t, std::size_t, std::size_t>
IntegrationGaussRegular<3>::getPosition(std::size_t nGauss, std::size_t igp)
{
    std::size_t gp_r = igp / (nGauss * nGauss);
    std::size_t gp_s = igp % (nGauss * nGauss);
    std::size_t gp_t = gp_s % nGauss;
    gp_s /= nGauss;
    return std::make_tuple(gp_r, gp_s, gp_t);
}

template <>
inline double IntegrationGaussRegular<1>::getPoint(std::size_t nGauss, std::size_t igp, double* x)
{
    auto pt = MathLib::GaussLegendre::getPoint(nGauss, igp);
    x[0] = pt.first;
    return pt.second;
}

template <>
inline double IntegrationGaussRegular<2>::getPoint(std::size_t nGauss, std::size_t igp, double* x)
{
    auto pos = getPosition(nGauss, igp);
    auto pt1 = MathLib::GaussLegendre::getPoint(nGauss, std::get<0>(pos));
    auto pt2 = MathLib::GaussLegendre::getPoint(nGauss, std::get<1>(pos));
    x[0] = pt1.first;
    x[1] = pt2.first;
    return pt1.second * pt2.second;
}

template <>
inline double IntegrationGaussRegular<3>::getPoint(std::size_t nGauss, std::size_t igp, double* x)
{
    auto pos = getPosition(nGauss, igp);
    auto pt1 = MathLib::GaussLegendre::getPoint(nGauss, std::get<0>(pos));
    auto pt2 = MathLib::GaussLegendre::getPoint(nGauss, std::get<1>(pos));
    auto pt3 = MathLib::GaussLegendre::getPoint(nGauss, std::get<2>(pos));
    x[0] = pt1.first;
    x[1] = pt2.first;
    x[2] = pt3.first;
    return pt1.second * pt2.second * pt3.second;
}

} //namespace

