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
#include <cstddef>

namespace GeoLib
{
/**
 *
 * \ingroup GeoLib
 *
 * \brief Class AABB is an axis aligned bounding box around a given
 * set of geometric points of (template) type PNT_TYPE.
 * */
template <typename PNT_TYPE = GeoLib::Point>
class AABB
{
public:
	/**
	 * construction of object, initialization the axis aligned bounding box
	 * */
	AABB() :
		_min_pnt(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()),
		_max_pnt(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min())
	{}

	/**
	 * copy constructor.
	 * @param src an axis aligned bounding box
	 */
	AABB(AABB<PNT_TYPE> const& src) :
		_min_pnt(src._min_pnt), _max_pnt(src._max_pnt)
	{}

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
			update(*it);
			it++;
		}
	}

	void update(PNT_TYPE const & pnt)
	{
		for (size_t k(0); k<3; k++) {
			if (pnt[k] < _min_pnt[k])
				_min_pnt[k] = pnt[k];
			if (_max_pnt[k] < pnt[k])
				_max_pnt[k] = pnt[k];
		}
	}

	/**
	 * check if point is in the axis aligned bounding box
	 * (employing containsPoint (double x, double y, double z))
	 */
	bool containsPoint(PNT_TYPE const & pnt) const
	{
		if (pnt[0] < _min_pnt[0] || _max_pnt[0] < pnt[0]) return false;
		if (pnt[1] < _min_pnt[1] || _max_pnt[1] < pnt[1]) return false;
		if (pnt[2] < _min_pnt[2] || _max_pnt[2] < pnt[2]) return false;
		return true;
	}

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
	 * @param other_aabb the AABB to test with
	 * @return true if the other AABB is contained in the AABB
	 * represented by this object
	 */
	bool containsAABB(AABB<PNT_TYPE> const& other_aabb) const
	{
		GeoLib::Point const& min_other(other_aabb.getMinPoint());
		GeoLib::Point const& max_other(other_aabb.getMaxPoint());
		for (unsigned k(0); k<3; k++) {
			if (_min_pnt[k] > min_other[k] || max_other[k] > _max_pnt[k])
				return false;
		}
		return true;
	}

protected:
	GeoLib::Point _min_pnt;
	GeoLib::Point _max_pnt;
private:
	void update(PNT_TYPE const * pnt)
	{
		update (*pnt);
	}
};
} // end namespace

#endif /* AXISALIGNEDBOUNDINGBOX_H_ */
