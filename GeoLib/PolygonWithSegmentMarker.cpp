/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PolygonWithSegmentMarker.h"

namespace GeoLib {
PolygonWithSegmentMarker::PolygonWithSegmentMarker(
    GeoLib::Polyline const& polyline)
    : GeoLib::Polygon(polyline, true),
      marker_(polyline.getNumberOfPoints(), false)
{
}

void PolygonWithSegmentMarker::markSegment(std::size_t seg_num, bool mark_val)
{
    marker_[seg_num] = mark_val;
}

bool PolygonWithSegmentMarker::isSegmentMarked(std::size_t seg_num) const
{
    return marker_[seg_num];
}

bool PolygonWithSegmentMarker::addPoint(std::size_t pnt_id)
{
    if (Polyline::addPoint(pnt_id)) {
        marker_.push_back(false);
        return true;
    }
    return false;
}

bool PolygonWithSegmentMarker::insertPoint(std::size_t pos, std::size_t pnt_id)
{
    if (Polyline::insertPoint(pos, pnt_id)) {
        marker_.insert(marker_.begin()+pos, marker_[pos]);
        return true;
    }
    return false;
}

}  // namespace GeoLib
