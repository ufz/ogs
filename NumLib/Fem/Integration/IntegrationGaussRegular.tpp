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
IntegrationGaussRegular<1>::getPosition(std::size_t /*nGauss*/, std::size_t igp)
{
    return {igp};
}

template <>
inline std::array<std::size_t, 2>
IntegrationGaussRegular<2>::getPosition(std::size_t nGauss, std::size_t igp)
{
    assert(igp < nGauss);
    return {igp / nGauss, igp % nGauss};
}

template <>
inline std::array<std::size_t, 3>
IntegrationGaussRegular<3>::getPosition(std::size_t nGauss, std::size_t igp)
{
    assert(igp < nGauss);
    std::size_t const gp_r = igp / (nGauss * nGauss);
    std::size_t const gp_s = igp % (nGauss * nGauss);
    return {gp_r, gp_s / nGauss, gp_s % nGauss };
}

template <>
inline MathLib::WeightedPoint1D IntegrationGaussRegular<1>::getWeightedPoint(std::size_t nGauss, std::size_t igp) //, double* x)
{
    assert(igp < nGauss);
	std::array<double,1> coords;
	coords[0] = MathLib::GaussLegendre::getPoint(nGauss, igp).first;
	return MathLib::WeightedPoint1D (coords,
			MathLib::GaussLegendre::getPoint(nGauss, igp).second);
}

template <>
inline MathLib::WeightedPoint2D IntegrationGaussRegular<2>::getWeightedPoint(std::size_t nGauss, std::size_t igp) //, double* x)
{
    assert(igp < nGauss);
    std::array<std::size_t, 2> const pos = getPosition(nGauss, igp);
    std::pair<double, double> const pt1 = MathLib::GaussLegendre::getPoint(nGauss, std::get<0>(pos));
    std::pair<double, double> const pt2 = MathLib::GaussLegendre::getPoint(nGauss, std::get<1>(pos));

    std::array<double, 2> coords;
    coords[0] = pt1.first;
    coords[1] = pt2.first;
    return MathLib::WeightedPoint2D (coords, pt1.second * pt2.second);
}

template <>
inline MathLib::WeightedPoint3D IntegrationGaussRegular<3>::getWeightedPoint(std::size_t nGauss, std::size_t igp) //, double* x)
{
    assert(igp < nGauss);
    std::array<std::size_t, 3> const pos = getPosition(nGauss, igp);

    std::pair<double, double> const pt1 = MathLib::GaussLegendre::getPoint(nGauss, std::get<0>(pos));
    std::pair<double, double> const pt2 = MathLib::GaussLegendre::getPoint(nGauss, std::get<1>(pos));
    std::pair<double, double> const pt3 = MathLib::GaussLegendre::getPoint(nGauss, std::get<2>(pos));

    std::array<double, 3> coords;
    coords[0] = pt1.first;
    coords[1] = pt2.first;
    coords[2] = pt3.first;

    return MathLib::WeightedPoint3D(coords, pt1.second * pt2.second * pt3.second);
}

} //namespace

