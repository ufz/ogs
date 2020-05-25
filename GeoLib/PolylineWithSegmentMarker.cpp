/**
 * \file
 * \author Thomas Fischer
 * \date   Apr 3, 2012
 * \brief  Implementation of the PolylineWithSegmentMarker class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PolylineWithSegmentMarker.h"

namespace GeoLib {
PolylineWithSegmentMarker::PolylineWithSegmentMarker(
    GeoLib::Polyline const& polyline)
    : GeoLib::Polyline(polyline), marker_(polyline.getNumberOfSegments(), false)
{
}

void PolylineWithSegmentMarker::markSegment(std::size_t seg_num, bool mark_val)
{
    marker_[seg_num] = mark_val;
}

bool PolylineWithSegmentMarker::isSegmentMarked(std::size_t seg_num) const
{
    return marker_[seg_num];
}

bool PolylineWithSegmentMarker::addPoint(std::size_t pnt_id)
{
    if (Polyline::addPoint(pnt_id)) {
        marker_.push_back(false);
        return true;
    }
    return false;
}

bool PolylineWithSegmentMarker::insertPoint(std::size_t pos, std::size_t pnt_id)
{
    if (Polyline::insertPoint(pos, pnt_id)) {
        marker_.insert(marker_.begin()+pos, marker_[pos]);
        return true;
    }
    return false;
}

}  // namespace GeoLib
