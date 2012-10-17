/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file AxisAlignedBoundingBox.h
 *
 * Created on 2010-04-22 by Thomas Fischer
 */

#ifndef AXISALIGNEDBOUNDINGBOX_H_
#define AXISALIGNEDBOUNDINGBOX_H_

#include "Point.h"
#include <limits>
#include <vector>

namespace GeoLib
{
/**
 *
 * \ingroup GeoLib
 *
 * \brief Class AABB is a bounding box around a given geometric entity
 * */
class AABB
{
public:
	/**
	 * construction of object, initialization the axis aligned bounding box
	 * */
	AABB ();

	/**
	 * copy constructor.
	 * @param src an axis aligned bounding box
	 * @return
	 */
	AABB(AABB const& src);

	/**
	 * construction of object using vector of points
	 * */
	AABB ( const std::vector<GeoLib::Point*>* points );

	template <typename Iterator>
	AABB(Iterator beg, Iterator end)
	{
		for (Iterator it(beg); it != end; it++) {
			for (std::size_t k(0); k<3; k++) {
				if ((*it)[k] < _min_pnt[k]) _min_pnt[k] = (*it)[k];
				if (_max_pnt[k] < (*it)[k]) _max_pnt[k] = (*it)[k];
			}
		}
	}

	void update (GeoLib::Point const & pnt);
	/**
	 * update axis aligned bounding box
	 */
	void update (double x, double y, double z);

	/**
	 * update axis aligned bounding box
	 */
	void update (const double* pnt)
	{
		update (pnt[0], pnt[1], pnt[2]);
	}

	/**
	 * check if point is in the axis aligned bounding box
	 * (employing containsPoint (double x, double y, double z))
	 */
	bool containsPoint (GeoLib::Point const & pnt,
	                    double eps = std::numeric_limits<double>::epsilon()) const;

	/**
	 * wrapper for GeoLib::Point
	 */
	bool containsPoint (const double* pnt, double eps =
	                            std::numeric_limits<double>::epsilon()) const;

	/**
	 * check if point described by its coordinates x, y, z is in
	 * the axis aligned bounding box
	 */
	bool containsPoint(double x, double y, double z, double eps =
					std::numeric_limits<double>::epsilon()) const;

	/**
	 * returns a point that coordinates are minimal for each dimension
	 * for the given point set
	 * @return a point
	 */
	GeoLib::Point const& getMinPoint () const { return _min_pnt; }

	/**
	 * returns a point that coordinates are maximal for each dimension
	 * within the given point set
	 * @return a point
	 */
	GeoLib::Point const& getMaxPoint () const { return _max_pnt; }

	/**
	 * Method checks if the given AABB object is contained within the
	 * AABB represented by this object.
	 * @param other the AABB to test with
	 * @return true if the other AABB is contained in the AABB
	 * represented by this object
	 */
	bool containsAABB (AABB const& other) const;

protected:
	GeoLib::Point _min_pnt;
	GeoLib::Point _max_pnt;
};
} // end namespace

#endif /* AXISALIGNEDBOUNDINGBOX_H_ */
