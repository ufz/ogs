// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    // clang-format off
    /**
     * \brief Write xdmf and h5 file with geometry and topology data.
     *
     * Nomenclature:
     *
     * <table>
     * <tr><th>OGS</th><th>XDMF</th></tr>
     * <tr><td></td><td></td></tr>
     * <tr><td>process variable</td><td>dynamic_attribute / variable attribute (legacy)</td></tr>
     * <tr><td>secondary variable, derived variable</td><td>dynamic attribute / variable attribute (legacy)</td></tr>
     * <tr><td>general time-dependent data in PropertyVector</td><td> dynamic attribute / variable attribute (legacy)</td></tr>
     * <tr><td>mesh nodes / geometry</td><td> static attribute</td></tr>
     * <tr><td>mesh elements topology</td><td> static attribute</td></tr>
     * <tr><td>MaterialIDs, bulk IDs, global IDs </td><td> static attribute</td></tr>
     * <tr><td>time-independent data in PropertyVector </td><td> static attribute</td></tr>
     * </table>
     *
     * @param meshes Meshes or NodePartitionedMeshes to be written to file(s)
     * @param filepath absolute or relative filepath to the hdf5 file
     * @param time_step number of the step (temporal collection)
     * @param initial_time time in seconds of the first time step
     * @param output_variable_names names of `PropertyVector`s that should be
     * outputted
     * @param static_attribute_names subset of output_variable_name that should
     * @anchor XdmfHdfWriter_static_attribute_names
     * be written as static attributes (written once, without a
     * time series) in addition to the always-static built-in names:
     * <ul><li>MaterialIDs</li><li>bulk_node_ids</li><li>bulk_element_ids</li>
     * <li>bulk_edge_ids</li>
     * <li>bulk_face_ids</li>
     * <li>global_node_ids</li>
     * <li>global_element_ids</li>
     * </ul>
     * @param use_compression if true, zlib compression in HDFWriter component
     * is used
     * @param n_files number of hdf5 output files
     * @param chunk_size_bytes Data will be split into chunks. The parameter
     * specifies the size (in bytes) of the largest chunk.
     * @param store_static_data_separately If true, geometry, topology, and
     * constant attributes are written to a separate static HDF5 file
     * (mesh_static.h5). If false, the single-file layout is used.
     */
    // clang-format on
    XdmfHdfWriter(
        std::vector<std::reference_wrapper<const MeshLib::Mesh>> const& meshes,
        std::filesystem::path const& filepath, unsigned long long time_step,
        double initial_time, std::set<std::string> const& output_variable_names,
        std::set<std::string> const& static_attribute_names,
        bool use_compression, unsigned int n_files,
        unsigned int chunk_size_bytes, bool store_static_data_separately);

    /**
     * \brief Adds data for either lazy (xdmf) or eager (hdf) writing algorithm
     * @param time time value of the current time_step
     */
    void writeStep(double time);

private:
    std::unique_ptr<HdfWriter> _static_hdf_writer;
    // hdf_writer must be destructed before xdmf_writer
    std::unique_ptr<HdfWriter> _hdf_writer;
    std::vector<std::unique_ptr<XdmfWriter>> _xdmf_writer;
};
}  // namespace MeshLib::IO
