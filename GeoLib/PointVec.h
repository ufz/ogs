/**
 * \file
 * \author Thomas Fischer / Karsten Rink
 * \date   2010-02-02
 * \brief  Definition of the PointVec class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// GeoLib
#include "AABB.h"
#include "Point.h"
#include "Station.h"

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
	 * @param points Pointer to a vector of pointers to GeoLib::Points.
	 * @attention{PointVec will take the ownership of (the pointer to)
	 * the vector, i.e. it deletes the vector in the destructor! The class
	 * takes also the ownership of the GeoLib::Points the pointers within
	 * the vector points at, i.e. it delete the points!}
	 * @param name_id_map A std::map that stores the relation name to point.
	 * @attention{PointVec will take the ownership of the vector, i.e. it
	 * deletes the names.}
	 * @param type the type of the point, \sa enum PointType
	 * @param rel_eps This is a relative error tolerance value for the test of identical points.
	 * The size of the axis aligned bounding box multiplied with the value of rel_eps gives the
	 * real tolerance \f$tol\f$. Two points \f$p_0, p_1 \f$ are identical iff
	 * \f$|p_1 - p_0| \le tol.\f$
	 */
	PointVec (const std::string& name, std::vector<Point*>* points,
	          std::map<std::string, std::size_t>* name_id_map = NULL,
	          PointType type = PointVec::POINT, double rel_eps = sqrt(std::numeric_limits<double>::epsilon()));

	/** Destructor deletes all Points of this PointVec. */
	virtual ~PointVec ();

	/**
	 * Method adds a Point to the (internal) standard vector and takes the ownership.
	 * If the given point is already included in the vector, the point will be destroyed and
	 * the id of the existing point will be returned.
	 * @param pnt the pointer to the Point
	 * @return the id of the point within the internal vector
	 */
	std::size_t push_back (Point* pnt);

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

	const std::vector<std::size_t>& getIDMap () const { return _pnt_id_map; }

	double getShortestPointDistance () const;
	const GeoLib::AABB<GeoLib::Point>& getAABB () const;

	/// Returns a subset of this point vector containing only the points specified in "subset" as PointWithID-objects
	std::vector<GeoLib::Point*>* getSubset(const std::vector<std::size_t> &subset);

private:
	/**
	 * Removes points out of the given point set that have the (nearly) same coordinates, i.e.
	 * the distance of two points is smaller than eps measured in the maximum norm.
	 * @param pnt_vec the given set of points stored in a vector
	 * @param pnt_id_map the id mapping
	 * @param eps if the distance (measured in maximum norm) between points \f$p_i\f$ and \f$p_j\f$
	 * is smaller than eps the points \f$p_i\f$ and \f$p_j\f$ are considered as equal.
	 */
	void makePntsUnique (std::vector<GeoLib::Point*>* pnt_vec, std::vector<std::size_t> &pnt_id_map, double eps = sqrt(std::numeric_limits<double>::min()));

	/**
	 * After the point set is modified (for example by makePntsUnique()) the mapping has to be corrected.
	 */
	void correctNameIDMapping();

	/** copy constructor doesn't have an implementation */
	// compiler does not create a (possible unwanted) copy constructor
	PointVec (const PointVec &);
	/** standard constructor doesn't have an implementation */
	// compiler does not create a (possible unwanted) standard constructor
	PointVec ();

	/** assignment operator doesn't have an implementation */
	// this way the compiler does not create a (possible unwanted) assignment operator
	PointVec& operator= (const PointVec& rhs);

	std::size_t uniqueInsert (Point* pnt);

	/** the type of the point (\sa enum PointType) */
	PointType _type;

	/**
	 * permutation of the geometric elements according
	 * to their lexicographical order
	 */
	std::vector<std::size_t> _pnt_id_map;

	/**
	 * method calculates the shortest distance of points inside the _pnt_vec
	 */
	void calculateShortestDistance ();
	/**
	 * squared shortest distance - calculated by calculateShortestAndLargestDistance, possible update by uniqueInsert
	 */
	double _sqr_shortest_dist;

	AABB<GeoLib::Point> _aabb;
};
} // end namespace

#endif /* POINTVEC_H_ */
