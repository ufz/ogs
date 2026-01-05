// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>
#include <vector>

namespace GeoLib
{
class GEOObjects;
}

namespace FileIO
{
namespace Legacy
{
/** Interface for handling geometry from OGS-5 and below (*.gli files) */

/** Reads geometric objects from file in gli format */
bool readGLIFileV4 (const std::string& fname,
                    GeoLib::GEOObjects& geo,
                    std::string& unique_name,
                    std::vector<std::string>& errors,
                    std::string const& gmsh_path);

/** Writes geometric objects from a specific geometry to a gli-file */
void writeGLIFileV4(const std::string& fname,
                    const std::string& geo_name,
                    const GeoLib::GEOObjects& geo);

/** Writes all geometric information to a gli-file */
void writeAllDataToGLIFileV4 (const std::string& fname, const GeoLib::GEOObjects& geo);

}  // namespace Legacy
} // end namespace FileIO
