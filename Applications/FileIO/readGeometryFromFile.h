// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>

namespace GeoLib
{
class GEOObjects;
}

namespace FileIO
{
void readGeometryFromFile(std::string const& fname,
                          GeoLib::GEOObjects& geo_objs,
                          std::string const& gmsh_path);
}
