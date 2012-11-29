/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file PolyLine.h
 *
 *  Created on 2010-01-14 by Thomas Fischer
 */

#ifndef POLYLINE_H_
#define POLYLINE_H_

// GeoLib
#include "GeoObject.h"
#include "Point.h"

// MathLib
#include "MathTools.h"

#include <cmath>
#include <vector>

namespace GeoLib
{
class Location
{
public:
	enum type {
		LEFT,
		RIGHT,
		BEYOND,
		BEHIND,
		BETWEEN,
		SOURCE,
		DESTINATION
	};
};

/**
 * \ingroup GeoLib
 *
 * \brief Class Polyline consists mainly of a reference to a point vector and
 * a vector that stores the indices in the point vector.
 * A polyline consists of at least one line segment. The polyline is specified by the points
 * of the line segments. The class Polyline stores ids of pointers to the points in the
 * _ply_pnt_ids vector.
 * */
class Polyline : public GeoObject
{
public:
	/** constructor
	 * \param pnt_vec a reference to the point vector
	 */
	Polyline(const std::vector<Point*>& pnt_vec);
	/**
	 * Copy constructor
	 * @param ply Polyline
	 */
	Polyline (const Polyline& ply);

	virtual ~Polyline() {}

	/** write the points to the stream */
	void write(std::ostream &os) const;

	/**
	 * Adds an id of a point at the end of the polyline. The id have to be inside
	 * the (internal) _ply_pnts vector the polyline is based on.
	 * @param pnt_id
	 */
	virtual void addPoint(std::size_t pnt_id);

	/**
	 * Method inserts a new point (that have to be inside the _ply_pnts vector)
	 * at the given position in the polyline.
	 * @param pos the position in the polyline, pos have to be a value into the interval [0, number of points)
	 * @param pnt_id the id of the new point in the vector of points the polyline is based on
	 */
	virtual void insertPoint(std::size_t pos, std::size_t pnt_id);

	/**
	 * Method removes a point from the polyline. The connecting line segments will
	 * be removed and the length of the polyline will be changed.
	 * @param pos a valid position within the polyline
	 */
	virtual void removePoint(std::size_t pos);

	/**
	 * Closes a polyline by adding a line segment that connects start- and end-point.
	 * \param ply A Polyline containing at least three points.
	 * \return A polygon.
	 */
	static Polyline* closePolyline(const Polyline& ply);

	/// Constructs one polyline from a vector of connected polylines.
	/// All polylines in this vector need to reference the same point vector.
	static Polyline* constructPolylineFromSegments(const std::vector<Polyline*> &ply_vec,
	                                               double prox = 0.0);

	/**
	 * returns the number of points,
	 * the number of segments is about one smaller
	 * */
	std::size_t getNumberOfPoints() const;

	/** returns true if the polyline is closed */
	bool isClosed() const;

	/**
	 * Method tests if the given id of a point (within the vector of points the polyline is
	 * based on) is inside the polyline.
	 * @param pnt_id the id of the point
	 * @return true if the point is part of the polyline, else false
	 */
	bool isPointIDInPolyline(std::size_t pnt_id) const;

	/**
	 * returns the index of the i-th polyline point
	 * in the point vector
	 */
	std::size_t getPointID(std::size_t i) const;

	/**
	 * Changes a point index for one point in a line
	 * @param idx Index of point in line
	 * @param id ID of point in PointVec object
	 */
	void setPointID(std::size_t idx, std::size_t id);

	/** \brief const access operator for the access to the i-th point of the polyline.
	 */
	const Point* operator[](std::size_t i) const;

	/**
	 * \brief returns the i-th point contained in the polyline
	 * */
	const Point* getPoint(std::size_t i) const;

	std::vector<Point*> const& getPointsVec () const;

	/**
	 * returns the length of the polyline until the k-th line segment
	 * @param k the k-th line segment
	 * @return the length of the polyline until the k-th line segment
	 */
	double getLength (std::size_t k) const;

	/**
	 * get the complete length vector
	 * @return the length vector of the polyline
	 */
	const std::vector<double>& getLengthVec () const;

	friend bool operator==(Polyline const& lhs, Polyline const& rhs);
protected:
	/**
	 * 2D method - ignores z coordinate. It calculates the location
	 * of the point relative to the k-th line segment of the polyline.
	 * (literature reference:
	 * Computational Geometry and Computer Graphics in C++; Michael J. Laszlo)
	 * @param k the number of line segment
	 * @param pnt the point
	 * @return a value of enum LOCATION
	 */
	Location::type getLocationOfPoint (std::size_t k, GeoLib::Point const & pnt) const;

	static bool pointsAreIdentical(const std::vector<Point*> &pnt_vec,
	                               std::size_t i,
	                               std::size_t j,
	                               double prox);

	/** a reference to the vector of pointers to the geometric points */
	const std::vector<Point*> &_ply_pnts;
	/** position of pointers to the geometric points */
	std::vector<std::size_t> _ply_pnt_ids;
	/**
	 * the k-th element of the vector contains the length of the polyline until the k-th segment
	 */
	std::vector<double> _length;
};

/** overload the output operator for class Polyline */
std::ostream& operator<< (std::ostream &os, Polyline const& pl);

bool containsEdge (const Polyline& ply, std::size_t id0, std::size_t id1);

bool isLineSegmentIntersecting (const Polyline& ply, GeoLib::Point const& s0, GeoLib::Point const& s1);

/**
 * comparison operator
 * @param lhs first polyline
 * @param rhs second polyline
 * @return true, if the polylines consists of the same sequence of line segments
 */
bool operator==(Polyline const& lhs, Polyline const& rhs);

} // end namespace

#endif /* POLYLINE_H_ */
