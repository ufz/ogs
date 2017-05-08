/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoMechanicsProcess-fwd.h"
#include "ThermoMechanicsProcess.h"

namespace ProcessLib
{
namespace ThermoMechanics
{

template class ThermoMechanicsProcess<2>;
template class ThermoMechanicsProcess<3>;

}   // namespace ThermoMechanics
}   // namespace ProcessLib
