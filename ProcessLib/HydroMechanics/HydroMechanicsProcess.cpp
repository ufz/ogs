/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HydroMechanicsProcess.h"

namespace ProcessLib
{
namespace HydroMechanics
{
template class HydroMechanicsProcess<2>;
template class HydroMechanicsProcess<3>;

}  // namespace HydroMechanics
}  // namespace ProcessLib
