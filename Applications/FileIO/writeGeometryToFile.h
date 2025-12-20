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

/// Write geometry given by the \c geo_objs object and specified by the name
/// stored in param \c geo_name either to a gml or a gli file. If the extension
/// given in the \c fname parameter is "gml" or "GML" a gml file is written. In
/// case the extension is "gli" or "GLI" a gli file is written.
void writeGeometryToFile(std::string const& geo_name,
    GeoLib::GEOObjects& geo_objs, std::string const& fname);
}
