/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PhaseFieldProcess-fwd.h"
#include "PhaseFieldProcess.h"

namespace ProcessLib
{
namespace PhaseField
{

template class PhaseFieldProcess<2>;
template class PhaseFieldProcess<3>;

}   // namespace PhaseField
}   // namespace ProcessLib
