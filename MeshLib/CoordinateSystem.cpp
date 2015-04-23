/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CoordinateSystem.h"

#include <limits>

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

namespace MeshLib
{

CoordinateSystem::CoordinateSystem(const Element &ele)
{
    GeoLib::AABB<MeshLib::Node> aabb(ele.getNodes(),
        ele.getNodes()+ele.getNNodes());
    MeshLib::CoordinateSystem coords(aabb);
    CoordinateSystem bboxCoordSys(getCoordinateSystem(aabb));
    if (bboxCoordSys.getDimension() >= ele.getDimension()) {
        _type = bboxCoordSys.getType();
    } else { // e.g. zero volume elements
        if (ele.getDimension()>=1)
            _type = CoordinateSystemType::X;
        if (ele.getDimension()>=2)
            _type |= CoordinateSystemType::Y;
        if (ele.getDimension()==3)
            _type |= CoordinateSystemType::Z;
    }
}

} // end
