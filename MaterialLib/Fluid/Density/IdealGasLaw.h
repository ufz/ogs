/*!
   \file  IdealGasLow.h
   \brief Declaration of class LIdealGasLow for fluid density by the ideal gas law
          depending on one variable linearly.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef IDEAL_GAS_LAW_H_
#define IDEAL_GAS_LAW_H_

#include <string>

#include "DensityType.h"

namespace MaterialLib
{

/// Fluid density by ideal gas low
class IdealGasLaw
{
    public:
        ///   \param molar_mass Molar mass of the gas phase.
        IdealGasLaw(const double molar_mass) : _molar_mass(molar_mass)
        {
        }

        /// Get density model name.
        std::string getName() const
        {
           return "Ideal gas law";  
        }

        DensityType getType() const
        {
            return DensityType::IDEAL_GAS;
        }

        /// Get density value
        /// \param T  Temperature in K.
        /// \param pg Gas phase pressure in Pa.        
        double getValue(const double T, const double pg) const
        {
           return _molar_mass * pg /(_gas_constant * T);  
        }
    private:       
       /// Normally denoted as R, unit \f$J(kmol K)^{-1}\f$.
       static constexpr double _gas_constant = 8315.41;
       
       /// Molar mass of gas phase.
       double _molar_mass;
};

} // end namespace
#endif

