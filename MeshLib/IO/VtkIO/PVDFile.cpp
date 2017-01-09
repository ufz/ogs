/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PVDFile.h"

#include <fstream>
#include <iomanip>
#include <limits>
#include <logog/include/logog.hpp>
#include "BaseLib/Error.h"

namespace MeshLib
{
namespace IO
{

void PVDFile::addVTUFile(const std::string &vtu_fname, double timestep)
{
    _datasets.push_back(std::make_pair(timestep, vtu_fname));

    std::ofstream fh(_pvd_filename.c_str());
    if (!fh) {
        OGS_FATAL("could not open file `%s'", _pvd_filename.c_str());
    }

    fh << std::setprecision(std::numeric_limits<double>::digits10);

    fh << "<?xml version=\"1.0\"?>\n"
           "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\""
           " compressor=\"vtkZLibDataCompressor\">\n"
           "  <Collection>\n";

    for (auto const& pair : _datasets)
        fh << "    <DataSet timestep=\"" << pair.first << "\" group=\"\" part=\"0\" file=\"" << pair.second << "\"/>\n";

    fh << "  </Collection>\n</VTKFile>\n";
}

} // IO
} // MeshLib
