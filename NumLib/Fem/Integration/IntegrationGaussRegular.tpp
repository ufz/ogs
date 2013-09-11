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

namespace detail
{

template <unsigned Dim, typename Method>
inline
MathLib::TemplateWeightedPoint<double, double, Dim>
aoeu(std::array<std::size_t, Dim> const& pos)
{
    std::array<double, Dim> coords;
    double weight = 1;
    for (unsigned d = 0; d < Dim; d++)
    {
        coords[d] = Method::X[pos[d]];
        weight *= Method::W[pos[d]];
    }

    return MathLib::TemplateWeightedPoint<double, double, Dim>(coords, weight);
}

}   // namespace detail

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

    switch (nGauss)
    {
        case 1: return detail::aoeu<1, MathLib::GaussLegendre<1>>(pos);
        case 2: return detail::aoeu<1, MathLib::GaussLegendre<2>>(pos);
        case 3: return detail::aoeu<1, MathLib::GaussLegendre<3>>(pos);
        case 4: return detail::aoeu<1, MathLib::GaussLegendre<4>>(pos);
    }

    return MathLib::WeightedPoint1D(std::array<double, 1>(), 0);
}

template <>
inline MathLib::WeightedPoint2D IntegrationGaussRegular<2>::getWeightedPoint(std::size_t nGauss, std::size_t igp)
{
    assert(igp < nGauss);
    std::array<std::size_t, 2> const pos = getPosition(nGauss, igp);

    switch (nGauss)
    {
        case 1: return detail::aoeu<2, MathLib::GaussLegendre<1>>(pos);
        case 2: return detail::aoeu<2, MathLib::GaussLegendre<2>>(pos);
        case 3: return detail::aoeu<2, MathLib::GaussLegendre<3>>(pos);
        case 4: return detail::aoeu<2, MathLib::GaussLegendre<4>>(pos);
    }

    return MathLib::WeightedPoint2D(std::array<double, 2>(), 0);
}

template <>
inline MathLib::WeightedPoint3D IntegrationGaussRegular<3>::getWeightedPoint(std::size_t nGauss, std::size_t igp)
{
    assert(igp < nGauss);
    std::array<std::size_t, 3> const pos = getPosition(nGauss, igp);

    switch (nGauss)
    {
        case 1: return detail::aoeu<3, MathLib::GaussLegendre<1>>(pos);
        case 2: return detail::aoeu<3, MathLib::GaussLegendre<2>>(pos);
        case 3: return detail::aoeu<3, MathLib::GaussLegendre<3>>(pos);
        case 4: return detail::aoeu<3, MathLib::GaussLegendre<4>>(pos);
    }

    return MathLib::WeightedPoint3D(std::array<double, 3>(), 0);
}

} //namespace

