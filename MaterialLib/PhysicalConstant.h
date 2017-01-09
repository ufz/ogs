/**
 * \copyright
 * Copyright (c) 2015-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

namespace MaterialLib
{
/**
 * Namespace containing physical constants
 * All members of this namespace should be given in SI standard units,
 * i.e., in terms of kg, m, s, K, mol, A, cd.
 */
namespace PhysicalConstant
{
/// Zero degrees Celsius in Kelvin
const double CelsiusZeroInKelvin = 273.15;

/**
  Ideal gas constant in SI standard units (J \per{mol} \per{K})

  - Source:               http://physics.nist.gov/cgi-bin/cuu/Value?r
  - Date visited:         2015-04-17
  - Standard uncertainty: 0.000 0075 J mol^-1 K^-1
*/
const double IdealGasConstant = 8.3144621;

/**
 * Molar masses of certain elements and chemical compounds
 */
namespace MolarMass
{
/**
 * Water.
 *
 * - Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
 * - Date visited: 2015-05-28
 *
 * According to the IUPAC report the molar mass of O is in the range [15.999 03,
 * 15.999 77] g/mol
 * and the molar mass of H is in the range [1.007 84, 1.008 11] g/mol
 */
const double Water = 0.018016;  ///< kg \per{mol}

/**
 * Nitrogen N<sub>2</sub>.
 *
 * - Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
 * - Date visited: 2015-05-28
 *
 * According to the IUPAC report the molar mass of N is in the range [14.006 43,
 * 14.007 28] g/mol
 */
const double N2 = 0.028013;  ///< kg \per{mol}

/**
 * Oxygen O<sub>2</sub>.
 *
 * - Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
 * - Date visited: 2015-05-28
 *
 * According to the IUPAC report the molar mass of O is in the range [15.999 03,
 * 15.999 77] g/mol
 */
const double O2 = 0.032;  ///< kg \per{mol}

/**
 * Air.
 *
 * - Source: http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
 */
const double Air = 0.02897;  ///< kg \per{mol}
}  // namespace MolarMass

/**
 * Specific gas constant \f$ R_\mathrm{s} = R/M \f$.
 * \$M\$ is the molar mass.
 * The unit of \f$ R_\mathrm{s} \f$ is J \per{kg} \per{K}.
 */
namespace SpecificGasConstant
{
/// Specific gas constant for water vapour.
const double WaterVapour = IdealGasConstant / MolarMass::Water;  // = 461.504;
}  // namespace SpecificGasConstant
}  // namespace PhysicalConstant
}  // namespace MaterialLib
