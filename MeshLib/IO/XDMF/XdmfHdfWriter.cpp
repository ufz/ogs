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

    std::vector<AttributeMeta> attributes = transformAttributes(mesh);
    std::vector<HdfData> hdf_data_attributes = {geometry.hdf, topology.hdf};
    std::transform(attributes.begin(), attributes.end(),
                   std::back_inserter(hdf_data_attributes),
                   [](AttributeMeta att) -> HdfData { return att.hdf; });

    std::function<bool(HdfData)> is_variable_hdf_attribute =
        [&variable_output_names](HdfData data) -> bool {
        return std::find(variable_output_names.begin(),
                         variable_output_names.end(),
                         data.name) != variable_output_names.end();
    };

    std::vector<HdfData> hdf_variable_attributes;
    std::copy_if(hdf_data_attributes.begin(), hdf_data_attributes.end(),
                 back_inserter(hdf_variable_attributes), is_variable_hdf_attribute);
    std::vector<HdfData> constant_attributes;
    std::copy_if(hdf_data_attributes.begin(), hdf_data_attributes.end(),
                 back_inserter(constant_attributes),
                 std::not_fn(is_variable_hdf_attribute));

    // HDF5
    std::filesystem::path const hdf_filepath =
        filepath.parent_path() / (filepath.stem().string() + ".h5");

    // geometry hdf data refers to the local data geometry and topology vector.
    // The hdf writer can write when data is out of scope.
    _hdf_writer = std::make_unique<HdfWriter>(std::move(constant_attributes),
                                              std::move(hdf_variable_attributes),
                                              time_step, hdf_filepath);
    // XDMF
    // The light data is only written by just one process
    if (!isFileManager())
    {
        return;
    }

    std::filesystem::path const xdmf_filepath =
        filepath.parent_path() / (filepath.stem().string() + ".xdmf");

    std::vector<XdmfData> xdmf_attributes;
    std::transform(attributes.begin(), attributes.end(),
                   std::back_inserter(xdmf_attributes),
                   [](AttributeMeta att) -> XdmfData { return att.xdmf; });

    std::function<bool(XdmfData)> is_variable_xdmf_attribute =
        [&variable_output_names](XdmfData data) -> bool {
        return std::find(variable_output_names.begin(),
                         variable_output_names.end(),
                         data.name) != variable_output_names.end();
    };

    std::vector<XdmfData> xdmf_variable_attributes;
    std::copy_if(xdmf_attributes.begin(), xdmf_attributes.end(),
                 back_inserter(xdmf_variable_attributes),
                 std::not1(is_variable_xdmf_attribute));

    std::vector<XdmfData> xdmf_constant_attributes;
    std::copy_if(xdmf_attributes.begin(), xdmf_attributes.end(),
                 back_inserter(xdmf_constant_attributes),
                 std::not_fn(is_variable_xdmf_attribute));

    _xdmf_writer = std::make_unique<Xdmf3Writer>(
        geometry.xdmf, topology.xdmf, std::move(xdmf_constant_attributes),
        std::move(xdmf_variable_attributes), xdmf_filepath, time_step);
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