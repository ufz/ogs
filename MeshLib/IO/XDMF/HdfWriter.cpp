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

static std::vector<Hdf5DimType> prependDimension(
    Hdf5DimType const prepend_value, std::vector<Hdf5DimType> const& dimensions)
{
    std::vector<Hdf5DimType> dims = {prepend_value};
    dims.insert(dims.end(), dimensions.begin(), dimensions.end());
    return dims;
}

static hid_t createDataSet(
    hid_t const data_type, std::vector<Hdf5DimType> const& data_dims,
    std::vector<Hdf5DimType> const& max_dims,
    [[maybe_unused]] std::vector<Hdf5DimType> const& chunk_dims,
    bool const use_compression, hid_t const section,
    std::string const& dataset_name)
{
    int const time_dim_local_size = data_dims.size() + 1;

    std::vector<Hdf5DimType> time_max_dims =
        prependDimension(H5S_UNLIMITED, max_dims);
    std::vector<Hdf5DimType> time_data_global_dims =
        prependDimension(1, max_dims);

    std::vector<Hdf5DimType> time_data_chunk_dims =
        prependDimension(1, chunk_dims);

    hid_t fspace =
        H5Screate_simple(time_dim_local_size, time_data_global_dims.data(),
                         time_max_dims.data());
    assert(fspace >= 0);

    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    assert(dcpl >= 0);

    hid_t status =
        H5Pset_chunk(dcpl, chunk_dims.size() + 1, time_data_chunk_dims.data());
    if (status < 0)
    {
        OGS_FATAL("H5Pset_layout failed for data set: {:s}.", dataset_name);
    }

    if (use_compression)
    {
        H5Pset_deflate(dcpl, default_compression_factor);
    }

    hid_t dataset = H5Dcreate2(section, dataset_name.c_str(), data_type, fspace,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);

    assert(dataset >= 0);
    H5Pclose(dcpl);
    assert(H5Sclose(fspace) >= 0);

    return dataset;
}
/**
 * \brief Assumes a dataset is already opened by createDatasetFunction
 * \details Defines what (nodes_data, data_type) will be written how (data
 * subsections: data_dims, offset_dims, max_dims, chunk_dims, time) where
 * (dataset and dataset_name)
 */
static void writeDataSet(
    void const* nodes_data,  // what
    hid_t const data_type,
    std::vector<Hdf5DimType> const& data_dims,  // how ...
    std::vector<Hdf5DimType> const& offset_dims,
    std::vector<Hdf5DimType> const& max_dims,
    [[maybe_unused]] std::vector<Hdf5DimType> const& chunk_dims,
    std::string const& dataset_name, int step, hid_t dataset)  // where
{
    Hdf5DimType hdf_step = step;
    Hdf5DimType time_steps = hdf_step + 1;

    std::vector<Hdf5DimType> time_data_local_dims = data_dims;
    std::vector<Hdf5DimType> time_max_dims =
        prependDimension(time_steps, max_dims);
    std::vector<Hdf5DimType> time_offsets =
        prependDimension(hdf_step, offset_dims);
    std::vector<hsize_t> count = prependDimension(1, time_data_local_dims);

    hid_t const io_transfer_property = createHDF5TransferPolicy();

    hid_t const mspace = H5Screate_simple(time_data_local_dims.size(),
                                          time_data_local_dims.data(), NULL);
    assert(H5Sselect_all(mspace) >= 0);

    hid_t status = H5Dset_extent(dataset, time_max_dims.data());
    if (status < 0)
    {
        OGS_FATAL("H5D set extent failed dataset '{:s}'.", dataset_name);
    }
    hid_t fspace = H5Dget_space(dataset);

    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, time_offsets.data(), NULL,
                        count.data(), NULL);

    status = H5Dwrite(dataset, data_type, mspace, fspace, io_transfer_property,
                      nodes_data);
    if (status < 0)
    {
        OGS_FATAL("H5Dwrite failed in dataset '{:s}'.", dataset_name);
    }

    H5Sclose(mspace);
    H5Pclose(io_transfer_property);

    return;
}
namespace MeshLib::IO
{
HdfWriter::HdfWriter(std::vector<HdfData> constant_attributes,
                     std::vector<HdfData>
                         variable_attributes,
                     int const initial_step,
                     std::filesystem::path const& filepath,
                     bool const use_compression)
    : _variable_attributes(std::move(variable_attributes)),
      _hdf5_filepath(filepath),
      _use_compression(checkCompression() && use_compression),
      _file(createFile(filepath)),
      _output_step(initial_step)
{
    _group = H5Gcreate2(_file, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    auto createAndWriteDataSet = [&](auto const& attribute) -> hid_t {
        hid_t dataset = createDataSet(
            attribute.data_type, attribute.data_space, attribute.file_space,
            attribute.chunk_space, _use_compression, _group, attribute.name);

        checkHdfStatus(dataset, "Creating HDF5 Dataset: {:s} failed.",
                       attribute.name);
        writeDataSet(attribute.data_start, attribute.data_type,
                     attribute.data_space, attribute.offsets,
                     attribute.file_space, attribute.chunk_space,
                     attribute.name, _output_step, dataset);
        return dataset;
    };

    for (auto const& attribute : constant_attributes)
    {
        hid_t dataset = createAndWriteDataSet(attribute);
        H5Dclose(dataset);
    }

    for (auto const& attribute : _variable_attributes)
    {
        hid_t dataset = createAndWriteDataSet(attribute);
        // datasets are kept open
        _datasets.insert({attribute.name, dataset});
    }
    _output_step++;
}

HdfWriter::~HdfWriter()
{
    for (auto& dataset : _datasets)
    {
        H5Dclose(dataset.second);
    }
    H5Gclose(_group);
    H5Fclose(_file);
}

void HdfWriter::writeStep()
{
    for (auto const& attribute : _variable_attributes)
    {
        auto const& dataset_hid = _datasets.find(attribute.name);
        if (dataset_hid == _datasets.end())
        {
            OGS_FATAL("Writing HDF5 Dataset: {:s} failed.", attribute.name);
        }

        writeDataSet(
            attribute.data_start, attribute.data_type, attribute.data_space,
            attribute.offsets, attribute.file_space, attribute.chunk_space,
            attribute.name, _output_step, _datasets.at(attribute.name));
    }
    _output_step++;
}
}  // namespace MeshLib::IO
