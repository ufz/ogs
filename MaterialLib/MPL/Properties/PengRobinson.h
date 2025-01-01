/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MaterialLib/PhysicalConstant.h"
#include "ParameterLib/Parameter.h"

namespace MaterialPropertyLib
{
/**
 * \brief Peng-Robinson equation of state
 *
 * This class implements the Peng-Robinson equation of state (PR-EOS),
 * a widely used cubic equation of state for describing the behaviour
 * of real gases, particularly hydrocarbons. It accounts for non-ideal
 * behaviour of fluids over a range of temperatures and pressures.
 *
 * \details
 * The equation is given in terms of the molar density \f$\rho\f$ as:
 * \f[
 * P = \frac{R T \rho}{1 - b \rho} - \frac{a \rho^2}{1 + 2b \rho - b^2 \rho^2}
 * \f]
 * where \f$P\f$ is the pressure, \f$T\f$ the temperature, \f$\rho\f$ the molar
 * density,
 * \f$R\f$ the universal gas constant, and \f$a\f$, \f$b\f$ are
 * substance-specific parameters.
 *
 * The parameters \f$a\f$ and \f$b\f$ are computed from the critical temperature
 * \f$T_c\f$ in Kelvin, critical pressure \f$p_c\f$ in Pascal, and
 * (dimensionless) acentric factor \f$\omega\f$ as follows: \f[ a = 0.457235
 * \frac{R^2 T_c^2}{p_c} \f] \f[ b = 0.077796 \frac{R T_c}{p_c} \f]
 *
 * The Peng-Robinson equation is applicable for a wide range of
 * substances (gases and liquids), particularly hydrocarbons, at conditions
 * ranging from subcritical to supercritical. The EOS is not suitable for solid
 * phases or very low-temperature applications where real gases behave ideally.
 *
 * All input parameters (temperature, pressure, density) are assumed to
 * be in SI units:
 * - Temperature \f$T\f$ in Kelvin [K]
 * - Pressure \f$P\f$ in Pascal [Pa]
 * - Mass density \f$\rho\f$ in \f$[kg/m^3]\f$
 * - The resulting properties will also follow SI units.
 *
 * Original source: D.-Y. Peng and D.B. Robinson, "A New Two-Constant
 * Equation of State," Industrial & Engineering Chemistry Fundamentals, vol. 15,
 * pp. 59-64, 1976.
 */
class PengRobinson final : public Property
{
public:
    explicit PengRobinson(const double Tc, const double pc, const double omega)
        : Tc_(Tc), pc_(pc), omega_(omega)
    {
        const double gas_constant =
            MaterialLib::PhysicalConstant::IdealGasConstant;
        a_ = 0.457235 * gas_constant * gas_constant * Tc_ * Tc_ / pc_;
        b_ = 0.077796 * gas_constant * Tc_ / pc_;
    }

    PropertyDataType value(
        MaterialPropertyLib::VariableArray const& variable_array,
        ParameterLib::SpatialPosition const& pos, double const t,
        double const dt) const override;
    PropertyDataType dValue(
        MaterialPropertyLib::VariableArray const& variable_array,
        Variable const variable, ParameterLib::SpatialPosition const& pos,
        double const t, double const dt) const override;

private:
    /// Critical temperature
    const double Tc_;
    /// Critical pressure
    const double pc_;
    /// Acentric factor
    const double omega_;
    /// Parameter a (cohesion pressure)
    double a_;
    /// Parameter b (covolume)
    double b_;
};
}  // namespace MaterialPropertyLib
