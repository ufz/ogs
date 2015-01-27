/*!
   \file  PressureDependentViscosity.h
   \brief Declaration of class for the pressure dependent viscosity model.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PRESSURE_DEPENDENT_VISCOSITY_H_
#define PRESSURE_DEPENDENT_VISCOSITY_H_

#include<string>
#include <algorithm>    // std::max

#include"ViscosityType.h"

namespace MaterialLib
{

/// Pressure dependent viscosity model,
class PressureDependentViscosity
{
    public:
        /*!
             \brief Viscosity defined by \f$\mu_0(1+\gamma(p-p_0))\f$
             \param mu0   Reference viscosity.
             \param p0    Reference pressure.
             \param gamma Paramater.
        */
        PressureDependentViscosity(const double mu0, const double p0,
                                   const double gamma)
                     :_mu0(mu0), _p0(p0), _gamma(gamma)
        {
        }

        /// Get model name.
        std::string getName() const
        {
           return "Pressure dependent viscosity" ;
        }

        ViscosityType getType() const
        {
            return ViscosityType::PRESSURE_DEPENDENT;
        }
        /// Get viscosity value
        /// \param p Pressure
        double getValue(const double p) const
        {
           return _mu0 * (1 + _gamma* (std::max(p, 0.) - _p0) );
        }
    private:
       double _mu0;    ///<  Reference viscosity..
       double _p0;     ///<  Reference pressure.
       double _gamma;  ///<  Paramater.
};

} // end namespace
#endif
