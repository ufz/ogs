/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
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
inline std::array<std::size_t, 1>
IntegrationGaussRegular<1>::getPosition(std::size_t nGauss, std::size_t igp)
{
    return {igp};
}

template <>
inline std::array<std::size_t, 2>
IntegrationGaussRegular<2>::getPosition(std::size_t nGauss, std::size_t igp)
{
    std::size_t gp_r = igp / nGauss;
    std::size_t gp_s = igp % nGauss;
    return {gp_r, gp_s};
}

template <>
inline std::array<std::size_t, 3>
IntegrationGaussRegular<3>::getPosition(std::size_t nGauss, std::size_t igp)
{
    std::size_t gp_r = igp / (nGauss * nGauss);
    std::size_t gp_s = igp % (nGauss * nGauss);
    std::size_t gp_t = gp_s % nGauss;
    gp_s /= nGauss;
    return {gp_r, gp_s, gp_t};
}

template <>
inline std::array<double, 2> IntegrationGaussRegular<1>::getWeightedPoint(std::size_t nGauss, std::size_t igp) //, double* x)
{
	return std::array<double, 2>({MathLib::GaussLegendre::getPoint(nGauss, igp).first,
		MathLib::GaussLegendre::getPoint(nGauss, igp).second});
//    auto pt = MathLib::GaussLegendre::getPoint(nGauss, igp);
//    x[0] = pt.first;
//    return pt.second;
}

template <>
inline std::array<double, 3> IntegrationGaussRegular<2>::getWeightedPoint(std::size_t nGauss, std::size_t igp) //, double* x)
{
    auto pos = getPosition(nGauss, igp);
    auto pt1 = MathLib::GaussLegendre::getPoint(nGauss, std::get<0>(pos));
    auto pt2 = MathLib::GaussLegendre::getPoint(nGauss, std::get<1>(pos));
//    x[0] = pt1.first;
//    x[1] = pt2.first;
    return std::array<double, 3>({pt1.first, pt2.first, pt1.second * pt2.second});
}

template <>
inline std::array<double, 4> IntegrationGaussRegular<3>::getWeightedPoint(std::size_t nGauss, std::size_t igp) //, double* x)
{
    auto pos = getPosition(nGauss, igp);
    auto pt1 = MathLib::GaussLegendre::getPoint(nGauss, std::get<0>(pos));
    auto pt2 = MathLib::GaussLegendre::getPoint(nGauss, std::get<1>(pos));
    auto pt3 = MathLib::GaussLegendre::getPoint(nGauss, std::get<2>(pos));
//    x[0] = pt1.first;
//    x[1] = pt2.first;
//    x[2] = pt3.first;
    return std::array<double, 4>({pt1.first, pt2.first, pt3.first, pt1.second * pt2.second * pt3.second});
}

} //namespace

