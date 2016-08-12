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
#include <algorithm>  // std::max

#include "BaseLib/ConfigTree.h"

#include "FluidDensityType.h"

namespace MaterialLib
{
namespace Fluid
{
/// Linear temperature dependent density model.
class LinearTemperatureDependentDensity
{
public:
    /// \param config  ConfigTree object which contains the input data
    ///                including <type>temperature_dependent</type> and it has
    ///                a tag of <liquid_density>
    LinearTemperatureDependentDensity(BaseLib::ConfigTree const* const config)
        :  //! \ogs_file_param{material__fluid__density__linear_temperature__rho0}
          _rho0(config->getConfigParameter<double>("rho0")),
          //! \ogs_file_param{material__fluid__density__linear_temperature__temperature0}
          _temperature0(config->getConfigParameter<double>("temperature0")),
          //! \ogs_file_param{material__fluid__density__linear_temperature__beta}
          _beta(config->getConfigParameter<double>("beta"))
    {
    }

    /// Get model name.
    std::string getName() const
    {
        return "Linear temperature dependent density";
    }

    FluidDensityType getType() const
    {
        return FluidDensityType::LINEAR_TEMPERATURE_DEPENDENT;
    }

    /// Get density value.
    /// \param T Temperature.
    double getValue(const double T) const
    {
        return _rho0 * (1 + _beta * (T - _temperature0));
    }

    /// Get the derivative of temperature.
    double getdValue(const double /* T */) const { return _rho0 * _beta; }
private:
    double _rho0;          ///<  Reference density.
    double _temperature0;  ///<  Reference temperature.
    double _beta;          ///<  Parameter.
};

}  // end namespace
}  // end namespace

#endif /* LINEARTEMPERATUREDEPENDENTDENSITY_H */
