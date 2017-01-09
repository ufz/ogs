/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CoordinateSystem.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

namespace MeshLib
{

CoordinateSystem::CoordinateSystem(const Element &ele)
{
    GeoLib::AABB const aabb(ele.getNodes(), ele.getNodes() + ele.getNumberOfNodes());
    CoordinateSystem const bboxCoordSys(getCoordinateSystem(aabb));
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

unsigned char CoordinateSystem::getCoordinateSystem(const GeoLib::AABB &bbox) const
{
    unsigned char coords = 0;

    const MathLib::Vector3 pt_diff(bbox.getMinPoint(), bbox.getMaxPoint());

    // The axis aligned bounding box is a from the right half-open interval.
    // Therefore, the difference between the particular coordinates of the
    // points is modified by the unit in the last place towards zero.
    if (std::nexttoward(std::abs(pt_diff[0]), 0.0) > .0)
        coords |= CoordinateSystemType::X;
    if (std::nexttoward(std::abs(pt_diff[1]), 0.0) > .0)
        coords |= CoordinateSystemType::Y;
    if (std::nexttoward(std::abs(pt_diff[2]), 0.0) > .0)
        coords |= CoordinateSystemType::Z;

    return coords;
}

} // end
