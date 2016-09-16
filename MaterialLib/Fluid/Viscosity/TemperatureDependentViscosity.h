/**
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   TemperatureDependentViscosity.h
 * \author: wenqing
 *
 *  Created on August 10, 2016, 10:45 AM
 */

#ifndef TEMPERATUREDEPENDENTVISCOSITY_H
#define TEMPERATUREDEPENDENTVISCOSITY_H

#include <string>
#include <cmath>

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/// Temperature dependent viscosity model,
class TemperatureDependentViscosity final : public FluidProperty
{
public:
    /// \param config  ConfigTree object which contains the input data
    ///                including <type>temperature_dependent</type> and it has
    ///                a tag of <viscosity>
    explicit TemperatureDependentViscosity(BaseLib::ConfigTree const& config)
        :  //! \ogs_file_param{material__fluid__viscosity__temperature_dependent__mu0}
    _mu0(config.getConfigParameter<double>("mu0")),
          //! \ogs_file_param{material__fluid__viscosity__temperature_dependent__tc}
    _temperature_c(config.getConfigParameter<double>("tc")),
          //! \ogs_file_param{material__fluid__viscosity__temperature_dependent__tv}
    _temperature_v(config.getConfigParameter<double>("tv"))
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Temperature dependent viscosity";
    }

    /// Get viscosity value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariableType.

    double getValue(const ArrayType& var_vals) const override
    {
        const double T = var_vals[static_cast<int> (PropertyVariableType::T)];
        return _mu0 * std::exp(-(T - _temperature_c) / _temperature_v);
    }

    /// Get the partial differential of the viscosity with respect to
    /// temperature.
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                   is given in enum class PropertyVariableType.
    /// \param var       Variable type.
    double getdValue(const ArrayType& var_vals,
            const PropertyVariableType var) const override
    {
        (void) var;
        const double T = var_vals[static_cast<int> (PropertyVariableType::T)];
        return -_mu0 * std::exp(-(T - _temperature_c) / _temperature_v);
    }

private:
    const double _mu0; ///<  Inital viscosity.
    const double _temperature_c; ///<  Reference temperature 1.
    const double _temperature_v; ///<  Reference temperature 2.
};

}  // end namespace
}  // end namespace

#endif /* TEMPERATUREDEPENDENTVISCOSITY_H */
