/**
 *  \brief Declare a class for dimensionless Gibbs free energy, region1.
 *
 *  \copyright
 *   Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   DimensionLessGibbsFreeEnergyRegion1.h
 *
 * Created on December 8, 2016, 12:31 PM
 */

#pragma once

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
class DimensionLessGibbsFreeEnergyRegion1
{
public:
    DimensionLessGibbsFreeEnergyRegion1() = default;

    /**
     * Get the value
     * @param pi  Dimension less temperature
     * @param tau Dimension less pressure
     * @return    The value
     */
    double get_gamma(const double tau, const double pi) const;

    /**
     * Get the 1st order partial derivative of the dimension less Gibbs free
     * energy with respect to dimension less temperature, tau
     *
     * @param pi  Dimension less temperature
     * @param tau Dimension less pressure
     * @return    The value
     */
    double get_dgamma_dtau(const double tau, const double pi) const;

    /**
     * Get the 2nd order partial derivative of the dimension less Gibbs free
     * energy with respect to dimension less temperature, tau
     *
     * @param pi  Dimension less temperature
     * @param tau Dimension less pressure
     * @return    The value
     */
    double get_dgamma_dtau_dtau(const double tau, const double pi) const;

    /**
     * Get the 1st order partial derivative of the dimension less Gibbs free
     * energy with respect to dimension less pressure, pi
     *
     * @param pi  Dimension less temperature
     * @param tau Dimension less pressure
     * @return    The value
     */
    double get_dgamma_dpi(const double tau, const double pi) const;

    /**
     * Get the 2nd order partial derivative of the dimension less Gibbs free
     * energy with respect to dimension less pressure, pi
     *
     * @param pi  Dimension less temperature
     * @param tau Dimension less pressure
     * @return    The value
     */
    double get_dgamma_dpi_dpi(const double tau, const double pi) const;

    /**
     * Get the 2nd order partial derivative of the dimension less Gibbs free
     * energy with respect to dimension less temperature, tau, and dimension
     * less pressure, pi
     * @param pi  Dimension less temperature
     * @param tau Dimension less pressure
     * @return    The value
     */
    double get_dgamma_dtau_dpi(const double tau, const double pi) const;
};

}  // end namespace
}  // end namespace
