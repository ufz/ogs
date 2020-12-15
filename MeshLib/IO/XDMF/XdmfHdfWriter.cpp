
#include "XdmfHdfWriter.h"

#include <algorithm>

#include "partition.h"
#include "transformData.h"

namespace MeshLib::IO
{
XdmfHdfWriter::XdmfHdfWriter(MeshLib::Mesh const& mesh,
                             std::filesystem::path const& filepath,
                             int const time_step)
{
    // transform Data into contiguous data and collect meta data
    Geometry const& geometry = transformGeometry(mesh);
    Topology const& topology = transformTopology(mesh, geometry.hdf.offsets[0]);
    std::vector<AttributeMeta> const attributes = transformAttributes(mesh);

    // HDF5
    std::filesystem::path const hdf_filepath =
        filepath.parent_path() / (filepath.stem().string() + ".h5");
    std::vector<HdfData> hdf_data_attributes;
    std::transform(attributes.begin(), attributes.end(),
                   std::back_inserter(hdf_data_attributes),
                   [](AttributeMeta att) -> HdfData { return att.hdf; });

    // geometry hdf data refers to the local data geometry and topology vector.
    // The hdf writer can write when data is out of scope.
    _hdf_writer = std::make_unique<HdfWriter>(geometry.hdf, topology.hdf,
                                              std::move(hdf_data_attributes),
                                              time_step, hdf_filepath);
    // XDMF
    // The light data is only written by just one process
    if (!isFileManager())
    {
        return;
    }

    std::filesystem::path const xdmf_filepath =
        filepath.parent_path() / (filepath.stem().string() + ".xdmf");
    std::vector<XdmfData> xdmf_data_attributes;
    std::transform(attributes.begin(), attributes.end(),
                   std::back_inserter(xdmf_data_attributes),
                   [](AttributeMeta att) -> XdmfData { return att.xdmf; });

    _xdmf_writer = std::make_unique<Xdmf3Writer>(
        geometry.xdmf, topology.xdmf, std::move(xdmf_data_attributes),
        xdmf_filepath, time_step);
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