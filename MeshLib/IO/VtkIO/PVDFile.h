/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

namespace MeshLib
{
namespace IO
{

/*! Writes a basic PVD file for use with Paraview.
 *
 */
class PVDFile
{
public:
    //! Set a PVD file path
    explicit PVDFile(std::string pvd_fname)
        : _pvd_filename(std::move(pvd_fname))
    {
    }

    //! Add a VTU file to this PVD file.
    void addVTUFile(std::string const& vtu_fname, double timestep);

private:
    std::string const _pvd_filename;
    std::vector<std::pair<double, std::string>> _datasets; // a vector of (time, VTU file name)
};

} // namespace IO
} // namespace MeshLib
