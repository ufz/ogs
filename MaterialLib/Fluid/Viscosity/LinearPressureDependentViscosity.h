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

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/// Linear pressure dependent viscosity model.
class LinearPressureDependentViscosity : public FluidProperty
{
public:
    /// \param config  ConfigTree object which contains the input data
    ///                including <type>linear_pressure</type> and it has
    ///                a tag of <viscosity>
    LinearPressureDependentViscosity(BaseLib::ConfigTree const& config)
        :  //! \ogs_file_param{material__fluid__viscosity__linear_pressure_dependent__mu0}
    _mu0(config.getConfigParameter<double>("mu0")),
          //! \ogs_file_param{material__fluid__viscosity__linear_pressure_dependent__p0}
    _p0(config.getConfigParameter<double>("p0")),
          //! \ogs_file_param{material__fluid__viscosity__linear_pressure_dependent__gamma}
    _gamma(config.getConfigParameter<double>("gamma"))
    {
    }

    /// Get model name.
    virtual std::string getName() const final
    {
        return "Linear pressure dependent viscosity";
    }

    /// Get viscosity value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariable.
    virtual double getValue(const double var_vals[]) const final
    {
        const double p = var_vals[static_cast<int> (PropertyVariable::pl)];
        return _mu0 * (1 + _gamma * (std::max(p,0.) - _p0));
    }

    /// Get the partial differential of the viscosity with respect to pressure.
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                  is given in enum class PropertyVariable.
    virtual double getdValue(const double /* var_vals */[],
                             const PropertyVariable /* var  */) const final
    {
        return _mu0 * _gamma;
    }

private:
    double _mu0;    ///<  Reference viscosity.
    double _p0;     ///<  Reference pressure.
    double _gamma;  ///<  Parameter.
};

}  // end namespace
}  // end namespace
#endif
