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

    std::vector<Hdf5DimType> const time_max_dims =
        prependDimension(H5S_UNLIMITED, max_dims);
    std::vector<Hdf5DimType> const time_data_global_dims =
        prependDimension(1, max_dims);

    std::vector<Hdf5DimType> const time_data_chunk_dims =
        prependDimension(1, chunk_dims);

    hid_t const fspace =
        H5Screate_simple(time_dim_local_size, time_data_global_dims.data(),
                         time_max_dims.data());
    assert(fspace >= 0);

    hid_t const dcpl = H5Pcreate(H5P_DATASET_CREATE);
    assert(dcpl >= 0);

    hid_t const status =
        H5Pset_chunk(dcpl, chunk_dims.size() + 1, time_data_chunk_dims.data());
    if (status < 0)
    {
        OGS_FATAL("H5Pset_layout failed for data set: {:s}.", dataset_name);
    }

    if (use_compression)
    {
        H5Pset_deflate(dcpl, default_compression_factor);
    }

    hid_t const dataset = H5Dcreate2(section, dataset_name.c_str(), data_type,
                                     fspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);

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
    std::string const& dataset_name, Hdf5DimType const step,
    hid_t const dataset)  // where
{
    Hdf5DimType const time_steps = step + 1;

    std::vector<Hdf5DimType> const time_data_local_dims = data_dims;
    std::vector<Hdf5DimType> const time_max_dims =
        prependDimension(time_steps, max_dims);
    std::vector<Hdf5DimType> const time_offsets =
        prependDimension(step, offset_dims);
    std::vector<hsize_t> const count =
        prependDimension(1, time_data_local_dims);

    hid_t const io_transfer_property = createHDF5TransferPolicy();

    hid_t const mspace = H5Screate_simple(time_data_local_dims.size(),
                                          time_data_local_dims.data(), NULL);
    assert(H5Sselect_all(mspace) >= 0);

    hid_t status = H5Dset_extent(dataset, time_max_dims.data());
    if (status < 0)
    {
        OGS_FATAL("H5D set extent failed dataset '{:s}'.", dataset_name);
    }
    hid_t const fspace = H5Dget_space(dataset);

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

/**
 * \brief Write vector with time values into open hdf file
 * \details In contrast to all other hdf write methods writing is only performed
 * by one process (is_file_manager_true). file handle is to an already opened
 * file
 */
static void writeTimeSeries(hid_t const file,
                            std::vector<double> const& step_times,
                            bool const is_file_manager)
{
    hsize_t const size = step_times.size();
    hid_t const memspace = H5Screate_simple(1, &size, NULL);
    hid_t const filespace = H5Screate_simple(1, &size, NULL);

    if (is_file_manager)
    {
        H5Sselect_all(memspace);
        H5Sselect_all(filespace);
    }
    else
    {
        H5Sselect_none(memspace);
        H5Sselect_none(filespace);
    }

    hid_t const dataset =
        H5Dcreate2(file, "/times", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT,
                   H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT,
             step_times.data());

    H5Dclose(dataset);
    H5Sclose(memspace);
    H5Sclose(filespace);
};
namespace MeshLib::IO
{
struct HdfWriter::HdfMesh final
{
    hid_t const group;
    std::string const name;
    std::map<std::string, hid_t> const datasets;
    std::vector<HdfData> const variable_attributes;
};

HdfWriter::HdfWriter(std::vector<MeshHdfData> meshes,
                     unsigned long long const initial_step,
                     std::filesystem::path const& filepath,
                     bool const use_compression,
                     bool const is_file_manager)
    : _hdf5_filepath(filepath),
      _file(createFile(filepath)),
      _meshes_group(
          H5Gcreate2(_file, "/meshes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)),
      _step_times{0},  // ToDo need to be initial time
      _use_compression(checkCompression() && use_compression),
      _is_file_manager(is_file_manager)
{
    for (auto const& mesh : meshes)
    {
        hid_t const group = H5Gcreate2(_meshes_group, mesh.name.c_str(),
                                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        auto const createAndWriteDataSet = [&](auto const& attribute) -> hid_t
        {
            hid_t const dataset = createDataSet(
                attribute.data_type, attribute.data_space, attribute.file_space,
                attribute.chunk_space, _use_compression, group, attribute.name);

            checkHdfStatus(dataset, "Creating HDF5 Dataset: {:s} failed.",
                           attribute.name);

            writeDataSet(attribute.data_start, attribute.data_type,
                         attribute.data_space, attribute.offsets,
                         attribute.file_space, attribute.chunk_space,
                         attribute.name, initial_step, dataset);
            return dataset;
        };

        for (auto const& attribute : mesh.constant_attributes)
        {
            hid_t const dataset = createAndWriteDataSet(attribute);
            H5Dclose(dataset);
        }

        std::map<std::string, hid_t> datasets;
        for (auto const& attribute : mesh.variable_attributes)
        {
            hid_t const dataset = createAndWriteDataSet(attribute);
            // datasets are kept open
            datasets.insert({attribute.name, dataset});
        }

        _hdf_meshes.push_back(std::make_unique<HdfMesh>(
            HdfMesh{group, mesh.name, datasets, mesh.variable_attributes}));
    }
}

HdfWriter::~HdfWriter()
{
    writeTimeSeries(_file, _step_times, _is_file_manager);

    for (auto const& mesh : _hdf_meshes)
    {
        for (auto const& dataset : mesh->datasets)
        {
            H5Dclose(dataset.second);
        }
        H5Gclose(mesh->group);
    }
    H5Gclose(_meshes_group);
    H5Fclose(_file);
}

void HdfWriter::writeStep(double const time)
{
    _step_times.push_back(time);
    auto const output_step = _step_times.size();

    for (auto const& mesh : _hdf_meshes)
    {
        for (auto const& attribute : mesh->variable_attributes)
        {
            auto const& dataset_hid = mesh->datasets.find(attribute.name);
            if (dataset_hid == mesh->datasets.end())
            {
                OGS_FATAL("Writing HDF5 Dataset: {:s} failed.", attribute.name);
            }

            writeDataSet(
                attribute.data_start, attribute.data_type, attribute.data_space,
                attribute.offsets, attribute.file_space, attribute.chunk_space,
                attribute.name, output_step, mesh->datasets.at(attribute.name));
        }
    }
}
}  // namespace MeshLib::IO
