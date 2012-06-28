/**
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

	GeoLib::Point const& getMinPoint () const { return _min_pnt; }
	GeoLib::Point const& getMaxPoint () const { return _max_pnt; }

protected:
	GeoLib::Point _min_pnt;
	GeoLib::Point _max_pnt;
};
} // end namespace

#endif /* AXISALIGNEDBOUNDINGBOX_H_ */
