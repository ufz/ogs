/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoHydroMechanicsProcess-fwd.h"
#include "ThermoHydroMechanicsProcess.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{

template class ThermoHydroMechanicsProcess<2>;
template class ThermoHydroMechanicsProcess<3>;

}   // namespace ThermoHydroMechanics
}   // namespace ProcessLib
