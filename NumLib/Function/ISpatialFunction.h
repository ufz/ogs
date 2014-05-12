/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef ISPATIALFUNCTION_H_
#define ISPATIALFUNCTION_H_

#include <vector>

#include "GeoLib/Point.h"

namespace NumLib
{

/**
 * \brief Interface class for any functions of spatial coordinates \f$f(x,y,z)\f$
 */
class ISpatialFunction
{
public:
	virtual ~ISpatialFunction(){}

	/**
	 * return a value at the given point
	 * \param pnt  a point object
	 * \return evaluated value
	 */
	virtual double operator()(const GeoLib::Point& pnt) const = 0;
};

} // NumLib

#endif //ISPATIALFUNCTION_H_
