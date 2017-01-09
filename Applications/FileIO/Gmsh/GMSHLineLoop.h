/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include <iosfwd>

namespace FileIO
{
namespace GMSH
{

class GMSHLine;

class GMSHLineLoop {
public:
    GMSHLineLoop(bool is_sfc=false);
    virtual ~GMSHLineLoop();
    void addLine(GMSHLine* line);
    bool isSurface() const { return _is_sfc; }
    void setSurface(bool is_sfc) { _is_sfc = is_sfc; }
    void write(std::ostream &os, std::size_t offset, std::size_t sfc_offset = 0) const;

private:
    std::vector<GMSHLine*> _lines;
    bool _is_sfc;
};

}  // end namespace GMSH
}  // end namespace FileIO
