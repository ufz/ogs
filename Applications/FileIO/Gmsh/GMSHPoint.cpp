/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GMSHPoint.h"

#include <cmath>
#include <limits>

namespace FileIO
{
namespace GMSH
{
GMSHPoint::GMSHPoint(GeoLib::Point const& pnt, std::size_t id,
                     double mesh_density)
    : GeoLib::Point(pnt, id), _mesh_density(mesh_density)
{
}

void GMSHPoint::write(std::ostream& os) const
{
    os << "Point(" << getID() << ") = {" << (*this)[0] << ", " << (*this)[1]
       << ", " << (*this)[2];
    if (fabs(_mesh_density) > std::numeric_limits<double>::epsilon())
    {
        os << ", " << _mesh_density << "};";
    }
    else
    {
        os << "};";
    }
}

std::ostream& operator<<(std::ostream& os, GMSHPoint const& p)
{
    p.write(os);
    return os;
}

}  // end namespace GMSH
}  // end namespace FileIO
