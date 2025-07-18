/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PVDFile.h"

#include <fstream>
#include <iomanip>
#include <limits>

#ifdef USE_PETSC
#include <mpi.h>

#include "BaseLib/MPI.h"
#endif

#include "BaseLib/Error.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"

namespace MeshLib
{
namespace IO
{
void PVDFile::addVTUFile(const std::string& vtu_fname, double timestep)
{
#ifdef USE_PETSC
    int mpi_size;
    MPI_Comm_size(BaseLib::MPI::OGS_COMM_WORLD, &mpi_size);
    if (mpi_size == 1)
    {
        _datasets.emplace_back(timestep, vtu_fname);
    }
    else
    {
        auto const vtu_file_name =
            getVtuFileNameForPetscOutputWithoutExtension(vtu_fname);

        _datasets.emplace_back(timestep, vtu_file_name + ".pvtu");
    }
#else
    _datasets.emplace_back(timestep, vtu_fname);
#endif

    std::ofstream fh(pvd_filename.c_str());
    if (!fh)
    {
        OGS_FATAL("could not open file `{:s}'", pvd_filename);
    }

    fh << std::setprecision(std::numeric_limits<double>::max_digits10);

    fh << "<?xml version=\"1.0\"?>\n"
          "<VTKFile type=\"Collection\" version=\"0.1\" "
          "byte_order=\"LittleEndian\""
          " compressor=\"vtkZLibDataCompressor\">\n"
          "  <Collection>\n";

    for (auto const& pair : _datasets)
    {
        fh << "    <DataSet timestep=\"" << pair.first
           << "\" group=\"\" part=\"0\" file=\"" << pair.second << "\"/>\n";
    }

    fh << "  </Collection>\n</VTKFile>\n";
}

}  // namespace IO
}  // namespace MeshLib
