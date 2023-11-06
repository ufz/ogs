/**
 * \file
 * \author Tobias Meisel
 * \date   2020-11-13
 * \brief  XdmfWriter which create contiguous data for geometry and topology
 * and writes this and all attributes to 1 xdmf + 1 hdf file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <set>

#include "HdfWriter.h"
#include "MeshLib/Mesh.h"
#include "XdmfWriter.h"

namespace MeshLib::IO
{
class XdmfHdfWriter final
{
public:
    /**
     * \brief Write xdmf and h5 file with geometry and topology data.
     * @param meshes Meshes or NodePartitionedMeshes to be written to file(s)
     * @param filepath absolute or relative filepath to the hdf5 file
     * @param time_step number of the step (temporal collection)
     * @param initial_time time in seconds of the first time step
     * @param variable_output_names names of all process variables (attributes)
     * that change over time
     * @param use_compression if true, zlib compression in HDFWriter component
     * is used
     * @param n_files number of hdf5 output files
     * @param chunk_size_bytes Data will be split into chunks. The parameter
     * specifies the size (in bytes) of the largest chunk.
     */
    XdmfHdfWriter(
        std::vector<std::reference_wrapper<const MeshLib::Mesh>> const& meshes,
        std::filesystem::path const& filepath, unsigned long long time_step,
        double initial_time, std::set<std::string> const& variable_output_names,
        bool use_compression, unsigned int n_files,
        unsigned int chunk_size_bytes);

    /**
     * \brief Adds data for either lazy (xdmf) or eager (hdf) writing algorithm
     * @param time time value of the current time_step
     */
    void writeStep(double time);

private:
    // hdf_writer must be destructed before xdmf_writer
    std::unique_ptr<HdfWriter> _hdf_writer;
    std::vector<std::unique_ptr<XdmfWriter>> _xdmf_writer;
};
}  // namespace MeshLib::IO
