/**
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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

class GMSHLineLoop final
{
public:
    explicit GMSHLineLoop(bool is_sfc = false);
    virtual ~GMSHLineLoop();
    bool isSurface() const { return is_sfc_; }
    void setSurface(bool is_sfc) { is_sfc_ = is_sfc; }
    void write(std::ostream& os, std::size_t offset,
               std::size_t sfc_offset = 0) const;

private:
    std::vector<GMSHLine*> lines_;
    bool is_sfc_;
};

}  // end namespace GMSH
}  // end namespace FileIO
