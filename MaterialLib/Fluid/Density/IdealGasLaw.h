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

#include "BaseLib/Error.h"

#include "MaterialLib/DensityType.h"
#include "MaterialLib/PhysicalConstant.h"

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
           return _molar_mass * pg /(PhysicalConstant::IdealGasConstant * T);  
        }

        /// Get the partial differential of density with the respect to
        /// or pressure.
        /// \param T  Temperature in K.
        /// \param pg Gas phase pressure in Pa.
        /// \param var_id Variable ID, 0 for temperature and 1 for pressure.        
        double getdValue(const double T, const double pg, const int var_id) const
        {
           if (var_id == 0)
               return -_molar_mass * pg /(PhysicalConstant::IdealGasConstant * T * T);
           else if (var_id == 1)
               return _molar_mass /(PhysicalConstant::IdealGasConstant * T);

           OGS_FATAL("Variable ID is larger than 1, "
                     "however the ideal gas law only has two variables." );
           return 0.;
        }

    private:
       /// Molar mass of gas phase.
       double _molar_mass;
};

} // end namespace
#endif

