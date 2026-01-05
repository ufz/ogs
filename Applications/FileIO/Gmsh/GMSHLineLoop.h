// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    ~GMSHLineLoop();
    bool isSurface() const { return _is_sfc; }
    void setSurface(bool is_sfc) { _is_sfc = is_sfc; }
    void write(std::ostream& os, std::size_t line_offset,
               std::size_t sfc_offset = 0) const;

private:
    std::vector<GMSHLine*> _lines;
    bool _is_sfc;
};

}  // end namespace GMSH
}  // end namespace FileIO
