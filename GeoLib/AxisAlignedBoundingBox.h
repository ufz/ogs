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
 * \brief Class AABB is an axis aligned bounding box around a given
 * set of geometric points.
 * */
class AABB
{
public:
	/**
	 * construction of object, initialization the axis aligned bounding box
	 * */
	AABB();

	/**
	 * copy constructor.
	 * @param src an axis aligned bounding box
	 */
	AABB(AABB const& src);

	/**
	 * Construction of object using input iterators. In contrast to give a vector
	 * this approach is more generic. You can use every (stl) container and
	 * C arrays as input for constructing the object. The constructor requires
	 * that the container stores pointers to point objects.
	 * @param first the input iterator to the initial position in the sequence
	 * @param last the input iterator to the final position in a sequence, i.e. [first, last)
	 */
	template <typename InputIterator>
	AABB(InputIterator first, InputIterator last)
	{
		InputIterator it(first);
		while (it != last) {
			this->update(*(reinterpret_cast<GeoLib::Point *const>(*it)));
			it++;
		}
	}

	void update(GeoLib::Point const & pnt);

	/**
	 * update axis aligned bounding box
	 */
	void update(const double* pnt)
	{
		update(pnt[0], pnt[1], pnt[2]);
	}

	/**
	 * check if point is in the axis aligned bounding box
	 * (employing containsPoint (double x, double y, double z))
	 */
	bool containsPoint(GeoLib::Point const & pnt,
	                   double eps = std::numeric_limits<double>::epsilon()) const;

	/**
	 * wrapper for GeoLib::Point
	 */
	bool containsPoint(const double* pnt, double eps =
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
	GeoLib::Point const& getMinPoint() const { return _min_pnt; }

	/**
	 * returns a point that coordinates are maximal for each dimension
	 * within the given point set
	 * @return a point
	 */
	GeoLib::Point const& getMaxPoint() const { return _max_pnt; }

	/**
	 * Method checks if the given AABB object is contained within the
	 * AABB represented by this object.
	 * @param other the AABB to test with
	 * @return true if the other AABB is contained in the AABB
	 * represented by this object
	 */
	bool containsAABB(AABB const& other) const;

protected:
	GeoLib::Point _min_pnt;
	GeoLib::Point _max_pnt;

private:
	/**
	 * update axis aligned bounding box
	 */
	void update(double x, double y, double z);
};
} // end namespace

#endif /* AXISALIGNEDBOUNDINGBOX_H_ */
