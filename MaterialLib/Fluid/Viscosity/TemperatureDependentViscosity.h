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

#include "ViscosityType.h"

namespace MaterialLib
{
namespace Fluid
{
/// Temperature dependent viscosity model,

class TemperatureDependentViscosity
{
public:
    /// \param config  ConfigTree object which contains the input data
    ///                including <type>temperature_dependent</type> and it has
    ///                a tag of <viscosity>
    TemperatureDependentViscosity(BaseLib::ConfigTree const* const config)
        :  //! \ogs_file_param{material__fluid__viscosity__temperature_dependent__mu0}
          _mu0(config->getConfigParameter<double>("mu0")),
          //! \ogs_file_param{material__fluid__viscosity__temperature_dependent__tc}
          _temperature_c(config->getConfigParameter<double>("tc")),
          //! \ogs_file_param{material__fluid__viscosity__temperature_dependent__tv}
          _temperature_v(config->getConfigParameter<double>("tv"))
    {
    }

    /// Get model name.
    std::string getName() const { return "Temperature dependent viscosity"; }
    ViscosityType getType() const
    {
        return ViscosityType::TEMPERATURE_DEPENDENT;
    }

    /// Get viscosity value
    /// \param T Temperature
    double getValue(const double T) const
    {
        return _mu0 * std::exp(-(T - _temperature_c) / _temperature_v);
    }

    /// Get the derivative of viscosity
    /// \param T Temperature
    double getdValue(const double T) const
    {
        return -_mu0 * std::exp(-(T - _temperature_c) / _temperature_v);
    }

private:
    double _mu0;            ///<  Inital viscosity.
    double _temperature_c;  ///<  Reference temperature 1.
    double _temperature_v;  ///<  Reference temperature 2.
};

}  // end namespace
}  // end namespace

#endif /* TEMPERATUREDEPENDENTVISCOSITY_H */
