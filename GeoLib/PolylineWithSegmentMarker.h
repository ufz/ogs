/*
 * PolylineWithSegmentMarker.h
 *
 *  Created on: Apr 3, 2012
 *      Author: fischeth
 */

#ifndef POLYLINEWITHSEGMENTMARKER_H_
#define POLYLINEWITHSEGMENTMARKER_H_

#include "Polyline.h"

namespace GeoLib {

class PolylineWithSegmentMarker: public GeoLib::Polyline {
public:
	PolylineWithSegmentMarker(GeoLib::Polyline const& polyline);
	virtual ~PolylineWithSegmentMarker();
	/**
	 * Method marks the segment (default mark is true).
	 * @param seg_num the segment number that should be marked
	 * @param mark_val the value of the flag (true or false)
	 */
	void markSegment(size_t seg_num, bool mark_val = true);
	/**
	 * Method returns the value of the mark for the given segment.
	 * @param seg_num segment number
	 * @return either true if the segment is marked or false else
	 */
	bool isSegmentMarked(size_t seg_num) const;

	/**
	 * Method calls @see Polyline::addPoint() and initializes the mark of the
	 * corresponding line segment.
	 * @param pnt_id @see Polyline::addPoint()
	 */
	virtual void addPoint(size_t pnt_id);

	/**
	 * Method calls the @see Polyline::insertPoint() and initializes the inserted line segment with the same
	 * value the previous line segment had.
	 * @param pos @see Polyline::insertPoint()
	 * @param pnt_id @see Polyline::insertPoint()
	 */
	virtual void insertPoint(size_t pos, size_t pnt_id);

private:
	std::vector<bool> _marker;
};

}

#endif /* POLYLINEWITHSEGMENTMARKER_H_ */
