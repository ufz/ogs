// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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

  Exact by definition since the 2019 revision of the SI: \f$R = N_A k_B\f$,
  and both the Avogadro constant and the Boltzmann constant are fixed
  exactly. There is therefore no standard uncertainty.

  - Source: SI Brochure, 9th edition (2019)
    https://www.bipm.org/en/publications/si-brochure
*/
constexpr double IdealGasConstant = 8.31446261815324;

/**
 * Boltzmann constant in SI standard units (J \per{K})
 *
 * Exact by definition since the 2019 revision of the SI: this value defines
 * the kelvin.
 *
 * - Source: SI Brochure, 9th edition (2019)
 *   https://www.bipm.org/en/publications/si-brochure
 */
constexpr double BoltzmannConstant = 1.380649e-23;  ///< J/K

/**
 * Planck constant in SI standard units (J s)
 *
 * Exact by definition since the 2019 revision of the SI: this value defines
 * the kilogram.
 *
 * - Source: SI Brochure, 9th edition (2019)
 *   https://www.bipm.org/en/publications/si-brochure
 */
constexpr double PlanckConstant = 6.62607015e-34;  ///< J s

/**
 * Critical point data of water.
 *
 * - Source: IAPWS R2-83(1992), "Release on the Values of Temperature, Pressure
 *   and Density of Ordinary and Heavy Water Substances at their Respective
 *   Critical Points".
 *   https://iapws.org/technical-guidance/release/crits
 */
namespace CriticalPoint
{
/// Critical temperature of ordinary water substance H<sub>2</sub>O.
/// Stated experimental uncertainty 0.1 K. The IAPWS-95 formulation adopts the
/// same value as the fixed reducing parameter of \f$\tau = T_c/T\f$, i.e. it is
/// not a quantity fitted by that formulation, see W. Wagner and A. Pruss, "The
/// IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water
/// Substance for General and Scientific Use", J. Phys. Chem. Ref. Data 31, 387
/// (2002), https://doi.org/10.1063/1.1461829.
constexpr double TemperatureWater = 647.096;  ///< K

/// Critical pressure of ordinary water substance H<sub>2</sub>O.
/// Stated experimental uncertainty 0.005 MPa. As for the critical
/// temperature, IAPWS-IF97 adopts this value as the upper bound of its
/// region 4 saturation line rather than reproducing it as a fitted quantity,
/// see IAPWS R7-97(2012), "Revised Release on the IAPWS Industrial
/// Formulation 1997 for the Thermodynamic Properties of Water and Steam",
/// Sec. 8.1. http://www.iapws.org/relguide/IF97-Rev.pdf
constexpr double PressureWater = 22.064e6;  ///< Pa
}  // namespace CriticalPoint

/**
 * Atomic masses of certain elements and molar masses of chemical compounds
 */
namespace MolarMass
{
/**
 * Water.
 *
 * Molar mass of ordinary water as adopted by the IAPWS-95 formulation. In OGS
 * this value only feeds \c SpecificGasConstant::WaterVapour, which is used
 * directly by \c WaterVapourDensity and, together with the IAPWS-IF97
 * Region-1 Gibbs derivative, by \c SaturationLuMcCartney::rho_w. The
 * IAPWS-IF97 property implementations themselves use the defined constant
 * \c SpecificGasConstant::WaterIF97 instead.
 *
 * - Source: IAPWS-95, Wagner and Pruss (2002), J. Phys. Chem. Ref. Data 31,
 *   387, Table 1.
 *   https://pubs.aip.org/aip/jpr/article/31/2/387/241937/The-IAPWS-Formulation-1995-for-the-Thermodynamic
 *   (the formulation tabulates this value). The primary metrology source
 *   and the uncertainty of this value are the IUPAC standard atomic weights,
 *   "Table of Standard Atomic Weights 2013" (as for \c N2 and \c O2 below),
 *   http://www.ciaaw.org/pubs/TSAW2013_xls.xls: its intervals for hydrogen
 *   [1.00784, 1.00811] and oxygen [15.99903, 15.99977] g/mol bound the
 *   molar mass of H<sub>2</sub>O to [18.01471, 18.01599] g/mol.
 */
constexpr double Water = 0.018015268;  ///< kg \per{mol}

/**
 * Nitrogen N<sub>2</sub>.
 *
 * - Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
 * - Date visited: 2015-05-28
 *
 * According to the IUPAC report the atomic mass of N is in the range
 * [14.00643, 14.00728] g/mol
 */
constexpr double N2 = 0.028013;  ///< kg \per{mol}

/**
 * Oxygen O<sub>2</sub>.
 *
 * - Source:       http://www.ciaaw.org/pubs/TSAW2013_xls.xls
 * - Date visited: 2015-05-28
 *
 * According to the IUPAC report the atomic mass of O is in the range
 * [15.99903, 15.99977] g/mol
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
/**
 * Specific gas constant for water vapour.
 *
 * Evaluates to 461.5231 J \per{kg} \per{K}. Note that this combines the exact
 * SI value of \c IdealGasConstant with the IAPWS-95 molar mass, and so differs
 * by 1.1e-5 in relative terms from the specific gas constant of IAPWS-95
 * itself (461.51805 J \per{kg} \per{K}, the value defined by IAPWS-95). That
 * difference is immaterial for the ideal-gas vapour relations using this
 * constant, but it must not be substituted for \c WaterIF97 below; note that
 * \c SaturationLuMcCartney::rho_w also consumes this constant together with
 * the IAPWS-IF97 Region-1 Gibbs derivative, so the difference does enter
 * results computed from that correlation.
 */
constexpr double WaterVapour = IdealGasConstant / MolarMass::Water;

/**
 * Specific gas constant of ordinary water as *defined* for IAPWS-IF97.
 *
 * This is not a measured physical constant and must not be replaced by
 * \c WaterVapour or updated to newer CODATA values: the regression
 * coefficients of IF97 were fitted with this value held fixed, so it is a
 * defined parameter of the formulation in the same sense as its reducing
 * quantities \f$p^*\f$ and \f$T^*\f$. Using any other value makes the IF97
 * equations disagree with their own published verification tables.
 *
 * - Source: IAPWS R7-97(2012), "Revised Release on the IAPWS Industrial
 *   Formulation 1997 for the Thermodynamic Properties of Water and Steam",
 *   Sec. 3 "Reference Constants", Eq. (1).
 *   http://www.iapws.org/relguide/IF97-Rev.pdf
 */
constexpr double WaterIF97 = 461.526;  ///< J \per{kg} \per{K}
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
