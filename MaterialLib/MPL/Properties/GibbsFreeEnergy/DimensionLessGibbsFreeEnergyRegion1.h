// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <boost/math/differentiation/autodiff.hpp>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

namespace MaterialLib
{
namespace Fluid
{
/**
 *  A class for dimensionless Gibbs free energy defined by
 *  \f[
 *     \gamma=\sum_{i=1}^{34}\left[
 *       n_i * (7.1-\pi)^{l_i}(\tau - 1.222)^{j_i}
 *     \right]
 *  \f]
 *  <a href="http://www.iapws.org/relguide/IF97-Rev.pdf">IF97-Rev</a>
 *
 *  Coefficients \f$n_i\f$, \f$j_i\f$ and \f$l_i\f$ are given in three static
 *  arrays in the cpp file.
 */
struct DimensionLessGibbsFreeEnergyRegion1
{
    /**
     * Get the value
     * @param pi  Dimension less temperature
     * @param tau Dimension less pressure
     * @return    The value
     */
    template <typename S, typename T>
    static boost::math::differentiation::promote<S, T> get_gamma(const S tau,
                                                                 const T pi);

    /**
     * Get the 1st order partial derivative of the dimension less Gibbs free
     * energy with respect to dimension less temperature, tau
     *
     * @param pi  Dimension less temperature
     * @param tau Dimension less pressure
     * @return    The value
     */
    template <typename S, typename T>
    static boost::math::differentiation::promote<S, T> get_dgamma_dtau(
        const S tau, const T pi);

    /**
     * Get the 2nd order partial derivative of the dimension less Gibbs free
     * energy with respect to dimension less temperature, tau
     *
     * @param pi  Dimension less temperature
     * @param tau Dimension less pressure
     * @return    The value
     */
    template <typename S, typename T>
    static boost::math::differentiation::promote<S, T> get_dgamma_dtau_dtau(
        const S tau, const T pi);

    /**
     * Get the 1st order partial derivative of the dimension less Gibbs free
     * energy with respect to dimension less pressure, pi
     *
     * @param pi  Dimension less temperature
     * @param tau Dimension less pressure
     * @return    The value
     */
    template <typename S, typename T>
    static boost::math::differentiation::promote<S, T> get_dgamma_dpi(
        const S tau, const T pi);

    /**
     * Get the 2nd order partial derivative of the dimension less Gibbs free
     * energy with respect to dimension less pressure, pi
     *
     * @param pi  Dimension less temperature
     * @param tau Dimension less pressure
     * @return    The value
     */
    template <typename S, typename T>
    static boost::math::differentiation::promote<S, T> get_dgamma_dpi_dpi(
        const S tau, const T pi);

    /**
     * Get the 2nd order partial derivative of the dimension less Gibbs free
     * energy with respect to dimension less temperature, tau, and dimension
     * less pressure, pi
     * @param pi  Dimension less temperature
     * @param tau Dimension less pressure
     * @return    The value
     */
    template <typename S, typename T>
    static boost::math::differentiation::promote<S, T> get_dgamma_dtau_dpi(
        const S tau, const T pi);
};

}  // namespace Fluid
}  // namespace MaterialLib
