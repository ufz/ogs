/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef LINESEGMENT_H_
#define LINESEGMENT_H_

#include <iostream>

#include "Point.h"

namespace GeoLib
{

struct LineSegment
{
	GeoLib::Point a;
	GeoLib::Point b;
};

std::ostream& operator<< (std::ostream& os, LineSegment const& s)
{
	os << "{(" << s.a << "), (" << s.b << ")}";
	return os;
}

std::ostream& operator<<(std::ostream& os,
                         std::pair<GeoLib::LineSegment const&,
                                   GeoLib::LineSegment const&> const& seg_pair)
{
	os << seg_pair.first << " x " << seg_pair.second;
	return os;
}

}
#endif
