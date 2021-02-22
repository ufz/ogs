/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "XdmfHdfWriter.h"

#include <algorithm>
#include <functional>

#include "partition.h"
#include "transformData.h"

namespace MeshLib::IO
{
XdmfHdfWriter::XdmfHdfWriter(MeshLib::Mesh const& mesh,
                             std::filesystem::path const& filepath,
                             int const time_step,
                             std::set<std::string>
                                 variable_output_names)
{
    // transform Data into contiguous data and collect meta data
    Geometry const& geometry = transformGeometry(mesh);
    Topology const& topology = transformTopology(mesh, geometry.hdf.offsets[0]);

    std::vector<AttributeMeta> all_attributes = transformAttributes(mesh);

    std::function<bool(AttributeMeta)> is_variable_attribute =
        [variable_output_names](AttributeMeta attributeMeta) {
            return std::find(variable_output_names.begin(),
                             variable_output_names.end(),
                             attributeMeta.hdf.name) !=
                   variable_output_names.end();
        };

    std::vector<AttributeMeta> variable_attributes;
    std::copy_if(all_attributes.begin(), all_attributes.end(),
                 back_inserter(variable_attributes), is_variable_attribute);
    std::vector<AttributeMeta> constant_attributes;
    std::copy_if(all_attributes.begin(), all_attributes.end(),
                 back_inserter(constant_attributes),
                 std::not1(is_variable_attribute));

    // HDF5
    std::filesystem::path const hdf_filepath =
        filepath.parent_path() / (filepath.stem().string() + ".h5");

    std::vector<HdfData> hdf_data_variable_attributes;
    std::transform(variable_attributes.begin(), variable_attributes.end(),
                   std::back_inserter(hdf_data_variable_attributes),
                   [](AttributeMeta att) -> HdfData { return att.hdf; });

    std::vector<HdfData> hdf_data_constant_attributes = {geometry.hdf,
                                                         topology.hdf};
    std::transform(constant_attributes.begin(), constant_attributes.end(),
                   std::back_inserter(hdf_data_constant_attributes),
                   [](AttributeMeta att) -> HdfData { return att.hdf; });

    // geometry hdf data refers to the local data geometry and topology vector.
    // The hdf writer can write when data is out of scope.
    _hdf_writer = std::make_unique<HdfWriter>(
        std::move(hdf_data_constant_attributes),
        std::move(hdf_data_variable_attributes), time_step, hdf_filepath);
    // XDMF
    // The light data is only written by just one process
    if (!isFileManager())
    {
        return;
    }

    std::filesystem::path const xdmf_filepath =
        filepath.parent_path() / (filepath.stem().string() + ".xdmf");

    std::vector<XdmfData> xdmf_data_variable_attributes;
    std::transform(variable_attributes.begin(), variable_attributes.end(),
                   std::back_inserter(xdmf_data_variable_attributes),
                   [](AttributeMeta att) -> XdmfData { return att.xdmf; });

    std::vector<XdmfData> xdmf_data_constant_attributes;
    std::transform(constant_attributes.begin(), constant_attributes.end(),
                   std::back_inserter(xdmf_data_constant_attributes),
                   [](AttributeMeta att) -> XdmfData { return att.xdmf; });

    _xdmf_writer = std::make_unique<Xdmf3Writer>(
        geometry.xdmf, topology.xdmf, std::move(xdmf_data_constant_attributes),
        std::move(xdmf_data_variable_attributes), xdmf_filepath, time_step);
}

void XdmfHdfWriter::writeStep(int const time_step, double const time) const
{
    _hdf_writer->writeStep(time_step);

    // XDMF
    // The light data is only written by just one process
    if (!_xdmf_writer)
    {
        return;
    }
    _xdmf_writer->writeStep(time_step, time);
}

}  // namespace MeshLib::IO