// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>
#include <utility>
#include <vector>

namespace MeshLib
{
namespace IO
{

/*! Writes a basic PVD file for use with Paraview.
 */
class PVDFile
{
public:
    //! Set a PVD file path
    explicit PVDFile(std::string pvd_fname)
        : pvd_filename(std::move(pvd_fname))
    {
    }

    //! Add a VTU file to this PVD file.
    void addVTUFile(std::string const& vtu_fname, double timestep);

    std::string const pvd_filename;

private:
    std::vector<std::pair<double, std::string>>
        _datasets;  // a vector of (time, VTU file name)
};

} // namespace IO
} // namespace MeshLib
