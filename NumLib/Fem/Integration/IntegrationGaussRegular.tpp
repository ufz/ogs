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
inline MathLib::WeightedPoint1D IntegrationGaussRegular<1>::getWeightedPoint(std::size_t nGauss, std::size_t igp)
{
    assert(igp < nGauss);
    std::array<std::size_t, 1> const pos = getPosition(nGauss, igp);

    std::array<double,1> coords;
    std::array<double,1> weight;
    switch (nGauss)
    {
        case 1:
            coords[0] = MathLib::GaussLegendre<1>::X[pos[0]];
            weight[0] = MathLib::GaussLegendre<1>::W[pos[0]];
            break;
        case 2:
            coords[0] = MathLib::GaussLegendre<2>::X[pos[0]];
            weight[0] = MathLib::GaussLegendre<2>::W[pos[0]];
            break;
        case 3:
            coords[0] = MathLib::GaussLegendre<3>::X[pos[0]];
            weight[0] = MathLib::GaussLegendre<3>::W[pos[0]];
            break;
        case 4:
            coords[0] = MathLib::GaussLegendre<4>::X[pos[0]];
            weight[0] = MathLib::GaussLegendre<4>::W[pos[0]];
            break;
        default:
            coords[0] = 0;
            weight[0] = 0;
    }

    return MathLib::WeightedPoint1D (coords, weight[0]);
}

template <>
inline MathLib::WeightedPoint2D IntegrationGaussRegular<2>::getWeightedPoint(std::size_t nGauss, std::size_t igp)
{
    assert(igp < nGauss);
    std::array<std::size_t, 2> const pos = getPosition(nGauss, igp);

    std::array<double, 2> coords;
    std::array<double, 2> weight;
    switch (nGauss)
    {
        case 1:
            coords[0] = MathLib::GaussLegendre<1>::X[pos[0]];
            coords[1] = MathLib::GaussLegendre<1>::X[pos[1]];
            weight[0] = MathLib::GaussLegendre<1>::W[pos[0]];
            weight[1] = MathLib::GaussLegendre<1>::W[pos[1]];
            break;
        case 2:
            coords[0] = MathLib::GaussLegendre<2>::X[pos[0]];
            coords[1] = MathLib::GaussLegendre<2>::X[pos[1]];
            weight[0] = MathLib::GaussLegendre<2>::W[pos[0]];
            weight[1] = MathLib::GaussLegendre<2>::W[pos[1]];
            break;
        case 3:
            coords[0] = MathLib::GaussLegendre<3>::X[pos[0]];
            coords[1] = MathLib::GaussLegendre<3>::X[pos[1]];
            weight[0] = MathLib::GaussLegendre<3>::W[pos[0]];
            weight[1] = MathLib::GaussLegendre<3>::W[pos[1]];
            break;
        case 4:
            coords[0] = MathLib::GaussLegendre<4>::X[pos[0]];
            coords[1] = MathLib::GaussLegendre<4>::X[pos[1]];
            weight[0] = MathLib::GaussLegendre<4>::W[pos[0]];
            weight[1] = MathLib::GaussLegendre<4>::W[pos[1]];
            break;
        default:
            coords[0] = 0;
            coords[1] = 0;
            weight[0] = 0;
            weight[1] = 0;
    }
    return MathLib::WeightedPoint2D (coords, weight[0]*weight[1]);
}

template <>
inline MathLib::WeightedPoint3D IntegrationGaussRegular<3>::getWeightedPoint(std::size_t nGauss, std::size_t igp)
{
    assert(igp < nGauss);
    std::array<std::size_t, 3> const pos = getPosition(nGauss, igp);

    std::array<double, 3> coords;
    std::array<double, 3> weight;
    switch (nGauss)
    {
        case 1:
            coords[0] = MathLib::GaussLegendre<1>::X[pos[0]];
            coords[1] = MathLib::GaussLegendre<1>::X[pos[1]];
            coords[2] = MathLib::GaussLegendre<1>::X[pos[2]];
            weight[0] = MathLib::GaussLegendre<1>::W[pos[0]];
            weight[1] = MathLib::GaussLegendre<1>::W[pos[1]];
            weight[2] = MathLib::GaussLegendre<1>::W[pos[2]];
            break;
        case 2:
            coords[0] = MathLib::GaussLegendre<2>::X[pos[0]];
            coords[1] = MathLib::GaussLegendre<2>::X[pos[1]];
            coords[2] = MathLib::GaussLegendre<2>::X[pos[2]];
            weight[0] = MathLib::GaussLegendre<2>::W[pos[0]];
            weight[1] = MathLib::GaussLegendre<2>::W[pos[1]];
            weight[2] = MathLib::GaussLegendre<2>::W[pos[2]];
            break;
        case 3:
            coords[0] = MathLib::GaussLegendre<3>::X[pos[0]];
            coords[1] = MathLib::GaussLegendre<3>::X[pos[1]];
            coords[2] = MathLib::GaussLegendre<3>::X[pos[2]];
            weight[0] = MathLib::GaussLegendre<3>::W[pos[0]];
            weight[1] = MathLib::GaussLegendre<3>::W[pos[1]];
            weight[2] = MathLib::GaussLegendre<3>::W[pos[2]];
            break;
        case 4:
            coords[0] = MathLib::GaussLegendre<4>::X[pos[0]];
            coords[1] = MathLib::GaussLegendre<4>::X[pos[1]];
            coords[2] = MathLib::GaussLegendre<4>::X[pos[2]];
            weight[0] = MathLib::GaussLegendre<4>::W[pos[0]];
            weight[1] = MathLib::GaussLegendre<4>::W[pos[1]];
            weight[2] = MathLib::GaussLegendre<4>::W[pos[2]];
            break;
        default:
            coords[0] = 0;
            coords[1] = 0;
            coords[2] = 0;
            weight[0] = 0;
            weight[1] = 0;
            weight[2] = 0;
    }

    return MathLib::WeightedPoint3D(coords, weight[0]*weight[1]*weight[2]);
}

} //namespace

