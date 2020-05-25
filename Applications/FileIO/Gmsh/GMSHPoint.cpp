/**
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cmath>
#include <limits>

#include "GMSHPoint.h"

namespace FileIO
{
namespace GMSH
{

GMSHPoint::GMSHPoint(GeoLib::Point const& pnt, std::size_t id, double mesh_density) :
    GeoLib::Point(pnt, id), mesh_density_(mesh_density)
{}

void GMSHPoint::write(std::ostream &os) const
{
    os << "Point(" << id_ << ") = {" << x_[0] << ", " << x_[1] << ", " << x_[2];
    if (fabs(mesh_density_) > std::numeric_limits<double>::epsilon()) {
        os << ", " << mesh_density_ << "};";
    } else {
        os << "};";
    }
}

std::ostream& operator<< (std::ostream &os, GMSHPoint const& p)
{
    p.write (os);
    return os;
}

}  // end namespace GMSH
}  // end namespace FileIO
