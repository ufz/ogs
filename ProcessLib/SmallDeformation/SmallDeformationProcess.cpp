/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SmallDeformationProcess.h"
#include "SmallDeformationProcess-impl.h"

namespace ProcessLib
{
namespace SmallDeformation
{

template class SmallDeformationProcess<2>;
template class SmallDeformationProcess<3>;

}   // namespace SmallDeformation
}   // namespace ProcessLib
