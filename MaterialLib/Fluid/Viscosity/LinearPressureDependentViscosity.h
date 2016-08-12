/*!
   \file  LinearPressureDependentViscosity.h
   \brief Declaration of class for the pressure dependent viscosity model.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
 */
#ifndef LINEAR_PRESSURE_DEPENDENT_VISCOSITY_H_
#define LINEAR_PRESSURE_DEPENDENT_VISCOSITY_H_

#include <string>
#include <algorithm>  // std::max

#include "BaseLib/ConfigTree.h"

#include "ViscosityType.h"

namespace MaterialLib
{
namespace Fluid
{
/// Linear pressure dependent viscosity model.
class LinearPressureDependentViscosity
{
public:
    /// \param config  ConfigTree object which contains the input data
    ///                including <type>linear_pressure</type> and it has
    ///                a tag of <viscosity>
    LinearPressureDependentViscosity(BaseLib::ConfigTree const* const config)
        :  //! \ogs_file_param{material__fluid__viscosity__linear_pressure_dependent__mu0}
          _mu0(config->getConfigParameter<double>("mu0")),
          //! \ogs_file_param{material__fluid__viscosity__linear_pressure_dependent__p0}
          _p0(config->getConfigParameter<double>("p0")),
          //! \ogs_file_param{material__fluid__viscosity__linear_pressure_dependent__gamma}
          _gamma(config->getConfigParameter<double>("gamma"))
    {
    }

    /// Get model name.
    std::string getName() const
    {
        return "Linear pressure dependent viscosity";
    }

    ViscosityType getType() const
    {
        return ViscosityType::LINEAR_PRESSURE_DEPENDENT;
    }

    /// Get viscosity value
    /// \param p Pressure
    double getValue(const double p) const
    {
        return _mu0 * (1 + _gamma * (std::max(p, 0.) - _p0));
    }

    /// Get the derivative of viscosity.
    /// \param p Pressure
    double getdValue(const double /* p */) const { return _mu0 * _gamma; }
private:
    double _mu0;    ///<  Reference viscosity.
    double _p0;     ///<  Reference pressure.
    double _gamma;  ///<  Parameter.
};

}  // end namespace
}  // end namespace
#endif
