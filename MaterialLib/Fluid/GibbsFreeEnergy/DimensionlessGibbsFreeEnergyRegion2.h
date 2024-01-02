/**
 *  \file
 *  \copyright
 *   Copyright (c) 2012-2024, OpenGeoSys Communty (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

namespace MaterialLib
{
namespace Fluid
{
/**
 *  The dimensionless Gibbs free energy in region2 defined in
 *  <a href="http://www.iapws.org/relguide/IF97-Rev.pdf">IF97-Rev</a>
 */
namespace DimensionlessGibbsFreeEnergyRegion2
{
/**
 * \param pi  Dimensionless temperature
 * \param tau Dimensionless pressure
 */
double getGamma(const double tau, const double pi);

/**
 * Get the 1st order partial derivative of the dimensionless Gibbs free
 * energy with respect to dimensionless temperature, tau
 */
double getdGammadTau(const double tau, const double pi);

/**
 * Get the 1st order partial derivative of the dimensionless Gibbs free
 * energy with respect to dimensionless pressure, pi
 */
double getdGammadPi(const double tau, const double pi);

}  // namespace DimensionlessGibbsFreeEnergyRegion2
}  // namespace Fluid
}  // namespace MaterialLib
