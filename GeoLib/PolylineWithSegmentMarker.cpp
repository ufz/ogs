/**
 * \file
 * \author Thomas Fischer
 * \date   Apr 3, 2012
 * \brief  Implementation of the PolylineWithSegmentMarker class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PolylineWithSegmentMarker.h"

namespace GeoLib {

PolylineWithSegmentMarker::PolylineWithSegmentMarker(GeoLib::Polyline const& polyline)
	: GeoLib::Polyline(polyline)
{
	const size_t n_pnts(getNumberOfPoints());
	_marker.resize(n_pnts);
	for (size_t k(0); k<n_pnts; k++) {
		_marker[k] = false;
	}
}

PolylineWithSegmentMarker::~PolylineWithSegmentMarker()
{}


void PolylineWithSegmentMarker::markSegment(size_t seg_num, bool mark_val)
{
	_marker[seg_num] = mark_val;
}
bool PolylineWithSegmentMarker::isSegmentMarked(size_t seg_num) const
{
	return _marker[seg_num];
}

void PolylineWithSegmentMarker::addPoint(size_t pnt_id)
{
	Polyline::addPoint(pnt_id);
	_marker.push_back(false);
}

void PolylineWithSegmentMarker::insertPoint(size_t pos, size_t pnt_id)
{
	Polyline::insertPoint(pos, pnt_id);
	_marker.insert(_marker.begin()+pos, _marker[pos]);
}

} // end GeoLib
