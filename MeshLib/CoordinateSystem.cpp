/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CoordinateSystem.h"

#include <cmath>

#include "GeoLib/AABB.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

namespace
{
unsigned char getCoordinateSystem(GeoLib::AABB const& bbox)
{
    unsigned char coords = 0;

    auto const [bbox_min, bbox_max] = bbox.getMinMaxPoints();
    Eigen::Vector3d const pt_diff = bbox_max - bbox_min;

    // The axis aligned bounding box is a from the right half-open interval.
    // Therefore, the difference between the particular coordinates of the
    // points is modified by the unit in the last place towards zero.
    if (std::nexttoward(std::abs(pt_diff[0]), 0.0) > .0)
    {
        coords |= MeshLib::CoordinateSystemType::X;
    }
    if (std::nexttoward(std::abs(pt_diff[1]), 0.0) > .0)
    {
        coords |= MeshLib::CoordinateSystemType::Y;
    }
    if (std::nexttoward(std::abs(pt_diff[2]), 0.0) > .0)
    {
        coords |= MeshLib::CoordinateSystemType::Z;
    }

    return coords;
}
}  // namespace

namespace MeshLib
{
CoordinateSystem::CoordinateSystem(const Element& ele)
{
    GeoLib::AABB const aabb(ele.getNodes(),
                            ele.getNodes() + ele.getNumberOfNodes());
    CoordinateSystem const bboxCoordSys(getCoordinateSystem(aabb));
    if (bboxCoordSys.getDimension() >= ele.getDimension())
    {
        _type = bboxCoordSys.getType();
    }
    else
    {  // e.g. zero volume elements
        if (ele.getDimension() >= 1)
        {
            _type = CoordinateSystemType::X;
        }
        if (ele.getDimension() >= 2)
        {
            _type |= CoordinateSystemType::Y;
        }
        if (ele.getDimension() == 3)
        {
            _type |= CoordinateSystemType::Z;
        }
    }
}

CoordinateSystem::CoordinateSystem(GeoLib::AABB const& bbox)
    : _type(getCoordinateSystem(bbox))
{
}

}  // namespace MeshLib
