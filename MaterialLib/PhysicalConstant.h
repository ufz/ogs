/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
constexpr double CelsiusZeroInKelvin = 273.15;

/**
  Ideal gas constant in SI standard units (J \per{mol} \per{K})

  - Source:               https://web.archive.org/web/20210603224431/https://physics.nist.gov/cgi-bin/cuu/Value?r
  - Date visited:         2015-04-17
  - Standard uncertainty: 0.000 0075 J mol^-1 K^-1
*/
constexpr double IdealGasConstant = 8.3144621;

/**
 * Atomic masses of certain elements and molar masses of chemical compounds
 */
namespace MolarMass
{
/**
 * Water.
 *
 * - Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
 * - Date visited: 2015-05-28
 *
 * According to the IUPAC report the atomic mass of O is in the range [15.999 03,
 * 15.999 77] g/mol
 * and the atomic mass of H is in the range [1.007 84, 1.008 11] g/mol
 */
constexpr double Water = 0.018016;  ///< kg \per{mol}

/**
 * Nitrogen N<sub>2</sub>.
 *
 * - Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
 * - Date visited: 2015-05-28
 *
 * According to the IUPAC report the atomic mass of N is in the range [14.006 43,
 * 14.007 28] g/mol
 */
constexpr double N2 = 0.028013;  ///< kg \per{mol}

/**
 * Oxygen O<sub>2</sub>.
 *
 * - Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
 * - Date visited: 2015-05-28
 *
 * According to the IUPAC report the atomic mass of O is in the range [15.999 03,
 * 15.999 77] g/mol
 */
constexpr double O2 = 0.032;  ///< kg \per{mol}

/**
 * Air.
 *
 * - Source: http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
 */
constexpr double Air = 0.02897;  ///< kg \per{mol}

/**
 * Hydrogen.
 *
 * - Source: https://pubchem.ncbi.nlm.nih.gov/compound/Hydrogen
 */
constexpr double H2 = 0.002016;  ///< kg \per{mol}
}  // namespace MolarMass

/**
 * Specific gas constant \f$ R_\mathrm{s} = R/M \f$.
 * \$M\$ is the molar mass.
 * The unit of \f$ R_\mathrm{s} \f$ is J \per{kg} \per{K}.
 */
namespace SpecificGasConstant
{
/// Specific gas constant for water vapour.
constexpr double WaterVapour = IdealGasConstant / MolarMass::Water;  // = 461.504;
}  // namespace SpecificGasConstant
/**
* Henry's law constant
* It is defined here as the ratio of the aqueous-phase concentration of a
* chemical to its equilibrium partial pressure in the gas phase.
*/
namespace HenryConstant
{
/**Henry constant for hydrogen, Please refer to
* De Nevers N. Physical and chemical equilibrium for chemical engineers[M].
* John Wiley & Sons, 2012.
*/
constexpr double HenryConstantH2 = 7.65e-6;  /// mol/Pa./m3
}  // namespace HenryConstant
}  // namespace PhysicalConstant
}  // namespace MaterialLib
