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

#include <cassert>

namespace NumLib
{

template <>
inline std::array<std::size_t, 1>
IntegrationGaussRegular<1>::getPositionIndices(std::size_t /*order*/, std::size_t igp)
{
    std::array<std::size_t, 1> result;
    result[0] = igp;
    return result;
}

template <>
inline std::array<std::size_t, 2>
IntegrationGaussRegular<2>::getPositionIndices(std::size_t order, std::size_t igp)
{
    assert(igp < order*order);
    std::array<std::size_t, 2> result;
    result[0] = igp / order;
    result[1] = igp % order;
    return result;
}

template <>
inline std::array<std::size_t, 3>
IntegrationGaussRegular<3>::getPositionIndices(std::size_t order, std::size_t igp)
{
    assert(igp < order*order*order);
    std::size_t const gp_r = igp / (order * order);
    std::size_t const gp_s = igp % (order * order);
    std::array<std::size_t, 3> result;
    result[0] = gp_r;
    result[1] = gp_s / order;
    result[2] = gp_s % order;
    return result;
}

template <std::size_t N_DIM>
inline
MathLib::TemplateWeightedPoint<double,double,N_DIM>
IntegrationGaussRegular<N_DIM>::getWeightedPoint(std::size_t order, std::size_t igp)
{
    assert(igp < std::pow(order, N_DIM));
    std::array<std::size_t, N_DIM> const pos = getPositionIndices(order, igp);

    switch (order)
    {
        case 1: return getWeightedPoint<MathLib::GaussLegendre<1>>(pos);
        case 2: return getWeightedPoint<MathLib::GaussLegendre<2>>(pos);
        case 3: return getWeightedPoint<MathLib::GaussLegendre<3>>(pos);
        case 4: return getWeightedPoint<MathLib::GaussLegendre<4>>(pos);
    }

    return MathLib::TemplateWeightedPoint<double, double, N_DIM>(std::array<double, N_DIM>(), 0);
}

template <std::size_t N_DIM>
template <typename Method>
inline
MathLib::TemplateWeightedPoint<double, double, N_DIM>
IntegrationGaussRegular<N_DIM>::getWeightedPoint(std::array<std::size_t, N_DIM> const& pos)
{
    std::array<double, N_DIM> coords;
    double weight = 1;
    for (unsigned d = 0; d < N_DIM; d++)
    {
        coords[d] = Method::X[pos[d]];
        weight *= Method::W[pos[d]];
    }

    return MathLib::TemplateWeightedPoint<double, double, N_DIM>(coords, weight);
}
} //namespace

