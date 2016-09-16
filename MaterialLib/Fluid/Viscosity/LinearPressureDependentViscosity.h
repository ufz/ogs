/*!
   \file  LinearPressureDependentViscosity.h
   \brief Declaration of class for the pressure dependent viscosity model.

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
 */
#ifndef LINEAR_PRESSURE_DEPENDENT_VISCOSITY_H_
#define LINEAR_PRESSURE_DEPENDENT_VISCOSITY_H_

#include <string>

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/// Linear pressure dependent viscosity model.
class LinearPressureDependentViscosity final : public FluidProperty
{
public:
    /// \param config  ConfigTree object which contains the input data
    ///                including <type>linear_pressure</type> and it has
    ///                a tag of <viscosity>
    explicit LinearPressureDependentViscosity(BaseLib::ConfigTree const& config)
        :  //! \ogs_file_param{material__fluid__viscosity__linear_pressure_dependent__mu0}
    _mu0(config.getConfigParameter<double>("mu0")),
          //! \ogs_file_param{material__fluid__viscosity__linear_pressure_dependent__p0}
    _p0(config.getConfigParameter<double>("p0")),
          //! \ogs_file_param{material__fluid__viscosity__linear_pressure_dependent__gamma}
    _gamma(config.getConfigParameter<double>("gamma"))
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Linear pressure dependent viscosity";
    }

    /// Get viscosity value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariableType.

    double getValue(const ArrayType& var_vals) const override
    {
        const double p = var_vals[static_cast<int> (PropertyVariableType::pl)];
        return _mu0 * (1 + _gamma * (std::max(p,0.) - _p0));
    }

    /// Get the partial differential of the viscosity with respect to pressure.
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                  is given in enum class PropertyVariableType.
    /// \param var       Variable type.
    double getdValue(const ArrayType& var_vals,
            const PropertyVariableType var) const override
    {
        (void) var_vals;
        (void) var;
        return _mu0 * _gamma;
    }

private:
    const double _mu0; ///<  Reference viscosity.
    const double _p0; ///<  Reference pressure.
    const double _gamma; ///<  Parameter.
};

}  // end namespace
}  // end namespace
#endif
