/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   LinearTemperatureDependentDensity.h
 * \author: wenqing
 *
 * Created on August 10, 2016, 11:34 AM
 */

#ifndef LINEARTEMPERATUREDEPENDENTDENSITY_H
#define LINEARTEMPERATUREDEPENDENTDENSITY_H

#include <string>

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/// Linear temperature dependent density model.
class LinearTemperatureDependentDensity final : public FluidProperty
{
public:
    /// \param config  ConfigTree object which contains the input data
    ///                including <type>temperature_dependent</type> and it has
    ///                a tag of <density>
    explicit LinearTemperatureDependentDensity(BaseLib::ConfigTree const& config)
        :  //! \ogs_file_param{material__fluid__density__linear_temperature__rho0}
    _rho0(config.getConfigParameter<double>("rho0")),
          //! \ogs_file_param{material__fluid__density__linear_temperature__temperature0}
    _temperature0(config.getConfigParameter<double>("temperature0")),
          //! \ogs_file_param{material__fluid__density__linear_temperature__beta}
    _beta(config.getConfigParameter<double>("beta"))
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Linear temperature dependent density";
    }

    /// Get density value.
    /// \param var_vals Variable values in an array. The order of its elements
    ///                 is given in enum class PropertyVariableType.
    double getValue(const ArrayType& var_vals) const override
    {
        const double T = var_vals[static_cast<int> (PropertyVariableType::T)];
        return _rho0 * (1 + _beta * (T - _temperature0));
    }

    /// Get the partial differential of the density with respect to temperature.
    /// \param var_vals  Variable values  in an array. The order of its elements
    ///                   is given in enum class PropertyVariableType.
    /// \param var       Variable type.
    double getdValue(const ArrayType& var_vals,
            const PropertyVariableType var) const override
    {
        (void) var_vals;
        (void) var;
        return _rho0 * _beta;
    }

private:
    const double _rho0; ///<  Reference density.
    const double _temperature0; ///<  Reference temperature.
    const double _beta; ///<  Parameter.
};

}  // end namespace
}  // end namespace

#endif /* LINEARTEMPERATUREDEPENDENTDENSITY_H */
