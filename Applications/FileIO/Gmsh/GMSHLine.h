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

#include <iosfwd>

namespace FileIO
{
namespace GMSH
{

class GMSHLine final {
public:
    GMSHLine(std::size_t start_point_id, std::size_t end_point_id);
    void write(std::ostream &os, std::size_t id) const;
    void resetLineData(std::size_t start_point_id, std::size_t end_point_id);

private:
    std::size_t _start_pnt_id;
    std::size_t _end_pnt_id;
};

} // end namespace GMSH
} // end namespace FileIO
