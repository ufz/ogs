/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com/LICENSE.txt
 *
 *
 * \file PointVec.h
 *
 * Created on 2010-02-02 by Thomas Fischer / Karsten Rink
 */

// GeoLib
#include "AxisAlignedBoundingBox.h"
#include "Point.h"
#include "Station.h"

// BaseLib
#include "binarySearch.h"
#include "quicksort.h"

#include <map>
#include <string>
#include <vector>

#ifndef POINTVEC_H_
#define POINTVEC_H_

#include "TemplateVec.h"

namespace GeoLib
{

class PointWithID;

/**
 * \ingroup GeoLib
 *
 * \brief This class manages pointers to Points in a std::vector along with a name.
 * It also handles the deleting of points. Additionally, each vector of points is identified by
 * a unique name from class GEOObject. For this reason PointVec should have
 * a name.
 * */
class PointVec : public TemplateVec<Point>
{
public:
	/// Signals if the vector contains object of type Point or Station
	enum PointType
	{
		POINT    = 0,
		STATION  = 1
	};

	/**
	 * Constructor initializes the name of the PointVec object,
	 * the internal pointer _pnt_vec to the raw points and the internal
	 * pointer the vector of names of the points
	 * and sets the type of PointVec.
	 * @param name the name of the point group
	 * @param points pointer to a vector of GeoLib::Pointers -
	 * PointVec will take the ownership of the vector,
	 * i.e. delete the points and the vector itself
	 * @param name_id_map the names to the points -
	 * PointVec will take the ownership of the vector, i.e. it
	 * deletes the names
	 * @param type the type of the point, \sa enum PointType
	 * @param rel_eps This is a relative error tolerance value for the test of identical points.
	 * The size of the axis aligned bounding box multiplied with the value of rel_eps gives the
	 * real tolerance \f$tol\f$. Two points \f$p_0, p_1 \f$ are identical iff
	 * \f$|p_1 - p_0| \le tol.\f$
	 * @return an object of type PointVec
	 */
	PointVec (const std::string& name, std::vector<Point*>* points,
	          std::map<std::string, size_t>* name_id_map = NULL,
	          PointType type = PointVec::POINT, double rel_eps = sqrt(std::numeric_limits<double>::min()));

	/** Destructor deletes all Points of this PointVec. */
	virtual ~PointVec ();

	/**
	 * Method adds a Point to the (internal) standard vector and takes the ownership.
	 * If the given point is already included in the vector, the point will be destroyed and
	 * the id of the existing point will be returned.
	 * @param pnt the pointer to the Point
	 * @return the id of the point within the internal vector
	 */
	size_t push_back (Point* pnt);

	/**
	 * push_back adds new elements at the end of the vector _data_vec.
	 * @param pnt a pointer to the point, PointVec takes ownership of the point
	 * @param name the name of the point
	 */
	virtual void push_back (Point* pnt, std::string const*const name);

	/**
	 * get the type of Point, this can be either POINT or STATION
	 *
	 */
	PointType getType() const { return _type; }

	std::vector<Point*>* filterStations(const std::vector<PropertyBounds> &bounds) const;

	const std::vector<size_t>& getIDMap () const { return _pnt_id_map; }

	double getShortestPointDistance () const;
	const GeoLib::AABB& getAxisAlignedBoundingBox () const;

	/// Creates a real copy of the point vector in memeory.
	static std::vector<GeoLib::Point*>* deepcopy(const std::vector<GeoLib::Point*> *pnt_vec);

	/// Returns a subset of this point vector containing only the points specified in "subset" as PointWithID-objects
	std::vector<GeoLib::Point*>* getSubset(const std::vector<size_t> &subset);

private:
	void makePntsUnique (std::vector<GeoLib::Point*>* pnt_vec, std::vector<size_t> &pnt_id_map, double eps = sqrt(std::numeric_limits<double>::min()));

	/** copy constructor doesn't have an implementation */
	// compiler does not create a (possible unwanted) copy constructor
	PointVec (const PointVec &);
	/** standard constructor doesn't have an implementation */
	// compiler does not create a (possible unwanted) standard constructor
	PointVec ();

	/** assignment operator doesn't have an implementation */
	// this way the compiler does not create a (possible unwanted) assignment operator
	PointVec& operator= (const PointVec& rhs);

	size_t uniqueInsert (Point* pnt);

	/** the type of the point (\sa enum PointType) */
	PointType _type;

	/**
	 * permutation of the geometric elements according
	 * to their lexicographical order
	 */
	std::vector<size_t> _pnt_id_map;

	/**
	 * method calculates the shortest distance of points inside the _pnt_vec
	 */
	void calculateShortestDistance ();
	/**
	 * squared shortest distance - calculated by calculateShortestAndLargestDistance, possible update by uniqueInsert
	 */
	double _sqr_shortest_dist;

	void calculateAxisAlignedBoundingBox ();
	AABB _aabb;
};
} // end namespace

#endif /* POINTVEC_H_ */
