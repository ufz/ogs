/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "HdfWriter.h"

#include <hdf5.h>

#include <string>
#include <utility>
#include <vector>

#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "fileIO.h"

template <typename... Args>
void checkHdfStatus(const hid_t status, std::string const& formatting,
                    Args&&... args)
{
    if (status < 0)
    {
        OGS_FATAL(formatting, std::forward<Args>(args)...);
    }
}

static unsigned short int const default_compression_factor = 1;

using namespace MeshLib::IO;

using namespace std::string_literals;

static std::string getTimeSection(int const step)
{
    return "t_"s + std::to_string(step);
}

static bool checkCompression()
{
    // Check if gzip compression is available and can be used for both
    // compression and decompression.
    if (htri_t avail = H5Zfilter_avail(H5Z_FILTER_DEFLATE); !avail)
    {
        WARN("gzip filter not available.\n");
        return false;
    }
    unsigned int filter_info;
    H5Zget_filter_info(H5Z_FILTER_DEFLATE, &filter_info);
    if (!(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ||
        !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED))
    {
        WARN("gzip filter not available for encoding and decoding.\n");
        return false;
    }
    return true;
}

static hid_t createStepGroup(hid_t const& file, int const step)
{
    std::string const& time_section = getTimeSection(step);

    // Open or create Group
    if (H5Lexists(file, time_section.c_str(), H5P_DEFAULT) > 0)
    {
        return H5Gopen2(file, time_section.c_str(), H5P_DEFAULT);
    }
    return H5Gcreate2(file, time_section.c_str(), H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
}

static hid_t writeDataSet(
    void const* nodes_data,  // what
    hid_t const data_type,
    std::vector<Hdf5DimType> const& data_dims,  // how ...
    std::vector<Hdf5DimType> const& dim_offsets,
    std::vector<Hdf5DimType> const& dim_maxs,
    [[maybe_unused]] std::vector<Hdf5DimType> const& chunk_dim,
    bool use_compression,
    hid_t const section,
    std::string const& dataset_name)  // where
{
    int const dim_size = data_dims.size();
    hid_t const memspace =
        H5Screate_simple(dim_size, data_dims.data(), nullptr);
    hid_t const filespace =
        H5Screate_simple(dim_size, dim_maxs.data(), nullptr);

    hid_t dataset_property = H5Pcreate(H5P_DATASET_CREATE);

    if (use_compression )
    {
        hid_t status =
            H5Pset_chunk(dataset_property, dim_size, dim_maxs.data());
        if (status != 0)
        {
            ERR("H5Pset_layout failed for data set: {:s}.", dataset_name);
        }
        H5Pset_deflate(dataset_property, default_compression_factor);
    }

    hid_t const dataset =
        H5Dcreate2(section, dataset_name.c_str(), data_type, filespace,
                   H5P_DEFAULT, dataset_property, H5P_DEFAULT);

    H5Pclose(dataset_property);
    H5Sclose(filespace);

    hid_t const dataset_filespace = H5Dget_space(dataset);

    std::vector<hsize_t> const stride(dim_size, 1);
    std::vector<hsize_t> const count(dim_size, 1);
    std::vector<hsize_t> const block = data_dims;

    hid_t status = H5Sselect_hyperslab(dataset_filespace, H5S_SELECT_SET,
                                 dim_offsets.data(), stride.data(),
                                 count.data(), block.data());
    if (status != 0)
    {
        ERR("H5Sselect_hyperslab failed in dataset '{:s}'.", dataset_name);
    }

    hid_t const io_transfer_property = createHDF5TransferPolicy();
    status = H5Dwrite(dataset, data_type, memspace, dataset_filespace,
                      io_transfer_property, nodes_data);
    if (status != 0)
    {
        ERR("H5Dwrite failed in dataset '{:s}'.", dataset_name);
    }

    H5Dclose(dataset);
    status = H5Sclose(memspace);
    H5Pclose(io_transfer_property);

    return (status >= 0 ? 1 : 0);
}
namespace MeshLib::IO
{
HdfWriter::HdfWriter(std::vector<HdfData> constant_attributes,
                     std::vector<HdfData>
                         variable_attributes,
                     int const step,
                     std::filesystem::path const& filepath,
                     bool const use_compression)
    : _variable_attributes(std::move(variable_attributes)),
      _hdf5_filepath(filepath),
      _use_compression(checkCompression() && use_compression),
      _file(createFile(filepath))
{

    std::string const& time_section = getTimeSection(step);
    hid_t const group_id = H5Gcreate2(_file, time_section.c_str(), H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);

    for (auto const& attribute : constant_attributes)
    {
        hid_t status = writeDataSet(attribute.data_start, attribute.data_type,
                                    attribute.data_space, attribute.offsets,
                                    attribute.file_space, attribute.chunk_space,
                                    _use_compression, group_id, attribute.name);

        checkHdfStatus(status, "Writing HDF5 Dataset: {:s} failed.",
                       attribute.name);
    }

    hid_t status = H5Gclose(group_id);

    checkHdfStatus(status, "HDF group could not be created!");
}

HdfWriter::~HdfWriter()
{
    H5Fclose(_file);
}

bool HdfWriter::writeStep(int const step) const
{
    hid_t const group = createStepGroup(_file, step);

    hid_t status = 0;
    for (auto const& attribute : _variable_attributes)
    {
        status = writeDataSet(attribute.data_start, attribute.data_type,
                              attribute.data_space, attribute.offsets,
                              attribute.file_space, attribute.chunk_space,
                              _use_compression, group, attribute.name);

        checkHdfStatus(status, "Writing HDF5 Dataset: {:s} failed.",
                       attribute.name);
    }

    H5Gclose(group);
    return (status >= 0);
}
}  // namespace MeshLib::IO
