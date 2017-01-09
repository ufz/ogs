/*!
   \file  IdealGasLaw.h
   \brief Declaration of class IdealGasLow for fluid density by the ideal gas
          law depending on one variable linearly.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
 */
#ifndef IDEAL_GAS_LAW_H_
#define IDEAL_GAS_LAW_H_

#include <cassert>
#include <string>

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialLib
{
namespace Fluid
{
/// Fluid density by ideal gas law
class IdealGasLaw final : public FluidProperty
{
public:
    ///   \param molar_mass Molar mass of the gas phase.
    explicit IdealGasLaw(const double molar_mass)
        : FluidProperty(), _molar_mass(molar_mass)
    {
    }

    /// Get density model name.
    std::string getName() const override { return "Ideal gas law"; }
    /// Get density value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariableType.
    double getValue(const ArrayType& var_vals) const override
    {
        return _molar_mass *
               var_vals[static_cast<int>(PropertyVariableType::p)] /
               (PhysicalConstant::IdealGasConstant *
                var_vals[static_cast<int>(PropertyVariableType::T)]);
    }

    /// Get the partial differential of the density with respect to temperature
    /// or gas pressure.
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                   is given in enum class PropertyVariableType.
    /// \param var       Variable type.
    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        const double p = var_vals[static_cast<int>(PropertyVariableType::p)];
        switch (var)
        {
            case PropertyVariableType::T:
                return dIdealGasLaw_dT(T, p);
            case PropertyVariableType::p:
                return dIdealGasLaw_dp(T, p);
            default:
                return 0.;
        }
    }

private:
    /// Molar mass of gas phase.
    const double _molar_mass;

    /// Get the partial differential of density with the respect to temperature
    /// \param T  Temperature in K.
    /// \param pg Gas phase pressure in Pa.
    double dIdealGasLaw_dT(const double T, const double pg) const
    {
        return -_molar_mass * pg / (PhysicalConstant::IdealGasConstant * T * T);
    }

    /// Get the partial differential of density with the respect to pressure
    /// \param T  Temperature in K.
    double dIdealGasLaw_dp(const double T, const double /* pg */) const
    {
        return _molar_mass / (PhysicalConstant::IdealGasConstant * T);
    }
};
}  // end namespace
}  // end namespace
#endif
