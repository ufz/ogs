/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef FILEIO_VTK_PVDFILE_H
#define FILEIO_VTK_PVDFILE_H

#include <string>
#include <iosfwd>

namespace FileIO
{

/*! Writes a basic PVD file for use with Paraview.
 *
 * The PVD file is automatically finished upon destruction
 * of this classes instance.
 */
class PVDFile
{
public:
    //! Create a new PVD file at the given path.
    explicit PVDFile(std::string const& fn);

    //! Add a VTU file to this PVD file.
    void addVTUFile(std::string const& fn, double timestep);

    //! Automatically finishes the PVD file's XML data.
    ~PVDFile();

private:
    std::ofstream _fh;
};

} // namespace FileIO

#endif // FILEIO_VTK_PVDFILE_H
