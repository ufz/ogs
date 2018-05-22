/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoMechanicalPhaseFieldProcess.h"
#include "ThermoMechanicalPhaseFieldProcess-impl.h"

namespace ProcessLib
{
namespace ThermoMechanicalPhaseField
{
template class ThermoMechanicalPhaseFieldProcess<2>;
template class ThermoMechanicalPhaseFieldProcess<3>;

}  // namespace ThermoMechanicalPhaseField
}  // namespace ProcessLib
