/*
 * \file PointVec.h
 *
 *  Created on: Feb 2, 2010
 *      Author: TF / KR
 */


// GEOLIB
#include "Point.h"
#include "Station.h"
#include "AxisAlignedBoundingBox.h"

// Base
#include "quicksort.h"
#include "binarySearch.h"

#include <vector>
#include <string>
#include <map>

#ifndef POINTVEC_H_
#define POINTVEC_H_

namespace GEOLIB {

/**
 * \ingroup GEOLIB
 *
 * \brief This class manages pointers to Points in a std::vector along with a name.
 * It also handles the deleting of points. Additionally, each vector of points is identified by
 * a unique name from class GEOObject. For this reason PointVec should have
 * a name.
 * */
class PointVec
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
	 * @param points pointer to a vector of GEOLIB::Pointers -
	 * PointVec will take the ownership of the vector,
	 * i.e. delete the points and the vector itself
	 * @param name_id_map the names to the points -
	 * PointVec will take the ownership of the vector, i.e. it
	 * deletes the names
	 * @param type the type of the point, \sa enum PointType
	 * @return an object of type PointVec
	 */
	PointVec (const std::string& name, std::vector<Point*>* points, std::map<std::string, size_t>* name_id_map = NULL,
				PointType type = PointVec::POINT);

	/** Destructor deletes all Points of this PointVec. */
	~PointVec ();

	/**
	 * Method adds a Point to the (internal) standard vector and takes the ownership.
	 * If the given point is already included in the vector, the point will be destroyed and
	 * the id of the existing point will be returned.
	 * @param pnt the pointer to the Point
	 * @return the id of the point within the internal vector
	 */
	size_t push_back (Point *pnt);

	/**
	 * push_back adds new elements at the end of the vector _pnt_vec.
	 * @param pnt a pointer to the point, PointVec takes ownership of the point
	 * @param name the name of the point
	 */
	void push_back (Point *pnt, const std::string& name);

	/**
	 * get the actual number of Points
	 */
	size_t size () const { return _pnt_vec->size(); }
	/**
	 * get the type of Point, this can be either POINT or STATION
	 *
	 */
	PointType getType() const { return _type; }

	/**
	 * getVector returns the internal vector of Points,
	 * you are not able to change the Points or the address of the vector.
	 */
	const std::vector<Point*>* getVector () const { return _pnt_vec; }

	std::vector<Point*> *filterStations(const std::vector<PropertyBounds> &bounds) const;

	/** sets the name of the object
	 * \param n the name as standard string */
	void setName(const std::string & n) { _name = n; }
	/** returns the name of the object */
	std::string getName () const { return _name; }

	/**
	 * search the vector of names for the ID of the point with the given name
	 * @param name the name of the point
	 * @param id the id of the point
	 * @return the id of the point
	 */
	bool getElementIDByName (const std::string& name, size_t &id) const;

	/**
	 * Method searchs for point with the given name. If it found a point with the name
	 * it returns a pointer to the point, else it returns the NULL pointer.
	 * @param name the name of the point
	 * @return the pointer to the point or NULL
	 */
	const Point* getElementByName (const std::string& name) const;

	/**
	 * The method returns true if the given element of type T
	 * can be found and the element has a name, else method returns false.
	 * @param data the data element, one wants to know the name
	 * @param name the name of the data element (if the data element is
	 * found and the data element has a name)
	 * @return if element is found and has a name: true, else false
	 */
	bool getNameOfElement (const Point* data, std::string& name) const;

	/**
	 * The method returns true if there is a name associated
	 * with the given id, else method returns false.
	 * @param id the id
	 * @param element_name if a name associated with the id
	 * is found name is assigned to element_name
	 * @return if there is name associated with the given id:
	 * true, else false
	 */
	bool getNameOfElementByID (size_t id, std::string& element_name) const;

	const std::vector<size_t>& getIDMap () const { return _pnt_id_map; }

	double getShortestPointDistance () const;
	const GEOLIB::AABB& getAxisAlignedBoundingBox () const;

private:
	void makePntsUnique (std::vector<GEOLIB::Point*>* pnt_vec, std::vector<size_t> &pnt_id_map);

	/** copy constructor doesn't have an implementation */
	// compiler does not create a (possible unwanted) copy constructor
	PointVec (const PointVec &);
	/** standard constructor doesn't have an implementation */
	// compiler does not create a (possible unwanted) standard constructor
	PointVec ();

	/** assignment operator doesn't have an implementation */
	// this way the compiler does not create a (possible unwanted) assignment operator
	PointVec& operator= (const PointVec& rhs);

	size_t uniqueInsert (Point *pnt);

	/**
	 * pointer to a vector of pointers to Points
	 *
	 * The destructor of PointVec will delete all GEOLIB::Points
	 * inside the vector.
	 */
	std::vector<Point*> *_pnt_vec;
	/**
	 * used to store the name associated with a point
	 */
	std::map<std::string, size_t>* _name_id_map;
	/** the type of the point (\sa enum PointType) */
	PointType _type;
	/** the name of the object */
	std::string _name;
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
	AABB aabb;
};

} // end namespace

#endif /* POINTVEC_H_ */
