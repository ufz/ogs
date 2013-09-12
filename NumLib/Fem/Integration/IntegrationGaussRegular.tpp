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

template <unsigned X, unsigned Y>
struct POW
{
    static unsigned const value = X * POW<X, Y>::value;
};

template <unsigned X>
struct POW<X, 0>
{
    static unsigned const value = 1;
};

template <unsigned Dim, unsigned Order, unsigned I, unsigned D>
struct position
{
    static unsigned const npoints = POW<Order, Dim - 1>::value;
    static unsigned const value =
        D == 0
        ? I / npoints
        : position<Dim - 1, Order, I % npoints, D - 1>::value;
};

template <unsigned Order, unsigned I>
struct position<1, Order, I, 0>
{
    static unsigned const npoints = Order;
    static unsigned const value = I;
};

template <>
inline std::array<std::size_t, 1>
IntegrationGaussRegular<1>::getPosition(std::size_t /*nGauss*/, std::size_t igp)
{
    return {{igp}};
}

template <>
inline std::array<std::size_t, 2>
IntegrationGaussRegular<2>::getPosition(std::size_t nGauss, std::size_t igp)
{
    assert(igp < nGauss*nGauss);
    return {{igp / nGauss, igp % nGauss}};
}

template <>
inline std::array<std::size_t, 3>
IntegrationGaussRegular<3>::getPosition(std::size_t nGauss, std::size_t igp)
{
    assert(igp < nGauss*nGauss*nGauss);
    std::size_t const gp_r = igp / (nGauss * nGauss);
    std::size_t const gp_s = igp % (nGauss * nGauss);
    return {{gp_r, gp_s / nGauss, gp_s % nGauss }};
}

template <std::size_t N_DIM>
inline
MathLib::TemplateWeightedPoint<double,double,N_DIM>
IntegrationGaussRegular<N_DIM>::getWeightedPoint(std::size_t nGauss, std::size_t igp)
{
    assert(igp < std::pow(nGauss,N_DIM));
    std::array<std::size_t, N_DIM> const pos = getPosition(nGauss, igp);

    switch (nGauss)
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

