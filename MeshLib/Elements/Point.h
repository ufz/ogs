/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHLIB_POINT_H_
#define MESHLIB_POINT_H_

#include "TemplateElement.h"
#include "PointRule1.h"

extern template class MeshLib::TemplateElement<MeshLib::PointRule1>;

namespace MeshLib
{
    using Point = TemplateElement<PointRule1>;
}

#endif  // MESHLIB_POINT_H_
