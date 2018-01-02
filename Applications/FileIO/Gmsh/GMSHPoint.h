/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

// GeoLib
#include "GeoLib/Point.h"

namespace FileIO
{
namespace GMSH
{
class GMSHPoint final : public GeoLib::Point
{
public:
    GMSHPoint(GeoLib::Point const& pnt, std::size_t id, double mesh_density);
    void write(std::ostream& os) const override;

private:
    double _mesh_density;
};

/** overload the output operator for class GMSHPoint */
std::ostream& operator<< (std::ostream &os, GMSHPoint const& p);

}  // end namespace GMSH
}  // end namespace FileIO
