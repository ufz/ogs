/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on August 16, 2019, 3:40 PM
 */

#pragma once

#include "MaterialLib/MPL/VariableType.h"  // for VariableArray
#include "ParameterLib/SpatialPosition.h"

namespace MaterialPropertyLib
{
class Phase;

/**
 * It gets the thermal expansion coefficient.
 *
 * If the the thermal expansion coefficient is given in the project file via
 * the media property of thermal_expansivity, e.g
 * \verbatim
 *     <property>
 *       <name>thermal_expansivity</name>
 *       <type>Constant</type>
 *       <value>2.07e-4</value>
 *     </property>
 * \endverbatim
 * it returns the value of the given property. Otherwise it returns the value
 * computed from the density model by the following formula
 * \f[
 *      (\frac{\partial \rho}{\partial T})/\rho
 * \f]
 * where \f$\rho\f$ is the density, \f$T\f$ is the temperature.
 */
double getLiquidThermalExpansivity(Phase const& phase,
                                   VariableArray const& vars,
                                   const double density,
                                   ParameterLib::SpatialPosition const& pos,
                                   double const t, double const dt);
}  // namespace MaterialPropertyLib
