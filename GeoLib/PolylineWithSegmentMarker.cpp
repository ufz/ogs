/**
 * \file
 * \author Thomas Fischer
 * \date   Apr 3, 2012
 * \brief  Implementation of the PolylineWithSegmentMarker class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PolylineWithSegmentMarker.h"

namespace GeoLib {
PolylineWithSegmentMarker::PolylineWithSegmentMarker(
    GeoLib::Polyline const& polyline)
    : GeoLib::Polyline(polyline), _marker(polyline.getNumberOfSegments(), false)
{
}

void PolylineWithSegmentMarker::markSegment(std::size_t seg_num, bool mark_val)
{
    _marker[seg_num] = mark_val;
}

bool PolylineWithSegmentMarker::isSegmentMarked(std::size_t seg_num) const
{
    return _marker[seg_num];
}

bool PolylineWithSegmentMarker::addPoint(std::size_t pnt_id)
{
    if (Polyline::addPoint(pnt_id)) {
        _marker.push_back(false);
        return true;
    }
    return false;
}

bool PolylineWithSegmentMarker::insertPoint(std::size_t pos, std::size_t pnt_id)
{
    if (Polyline::insertPoint(pos, pnt_id)) {
        _marker.insert(_marker.begin()+pos, _marker[pos]);
        return true;
    }
    return false;
}

} // end GeoLib
