/*!
   \file  IdealGasLow.h
   \brief Declaration of class LIdealGasLow for fluid density by the ideal gas
   law
          depending on one variable linearly.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
 */
#ifndef IDEAL_GAS_LAW_H_
#define IDEAL_GAS_LAW_H_

#include <string>

#include "FluidDensityType.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialLib
{
namespace Fluid
{
/// Fluid density by ideal gas low

class IdealGasLaw
{
public:
    ///   \param molar_mass Molar mass of the gas phase.
    IdealGasLaw(const double molar_mass)
        : _molar_mass(molar_mass),
          _derivative_functions{&IdealGasLaw::dIdealGasLaw_dT,
                                &IdealGasLaw::dIdealGasLaw_dp}
    {
    }

    /// Get density model name.

    std::string getName() const { return "Ideal gas law"; }
    FluidDensityType getType() const { return FluidDensityType::IDEAL_GAS; }
    /// Get density value
    /// \param T  Temperature in K.
    /// \param pg Gas phase pressure in Pa.

    double getValue(const double T, const double pg) const
    {
        return _molar_mass * pg / (PhysicalConstant::IdealGasConstant * T);
    }

    /// Get the partial differential of density with the respect to
    /// or pressure.
    /// \param T  Temperature in K.
    /// \param pg Gas phase pressure in Pa.
    /// \param var_id Variable ID, 0 for temperature and 1 for pressure.

    double getdValue(const double T, const double pg, const int var_id) const
    {
        assert(var_id > -1 && var_id < 2);

        return (this->*_derivative_functions[var_id])(T, pg);
    }

private:
    /// Molar mass of gas phase.
    double _molar_mass;

    /// Get the partial differential of density with the respect to temperature
    /// \param T  Temperature in K.
    /// \param pg Gas phase pressure in Pa.

    double dIdealGasLaw_dT(const double T, const double pg) const
    {
        return -_molar_mass * pg / (PhysicalConstant::IdealGasConstant * T * T);
    }

    /// Get the partial differential of density with the respect to pressure
    /// \param T  Temperature in K.
    /// \param pg Gas phase pressure in Pa.

    double dIdealGasLaw_dp(const double T, const double /* pg */) const
    {
        return _molar_mass / (PhysicalConstant::IdealGasConstant * T);
    }

    typedef double (IdealGasLaw::*ptr2_derivative_f)(const double,
                                                     const double) const;

    /// An array of pointer to derivative functions.
    ptr2_derivative_f _derivative_functions[2];
};
}  // end namespace
}  // end namespace
#endif
