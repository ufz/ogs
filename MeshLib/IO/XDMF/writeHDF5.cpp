/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "writeHDF5.h"

#include <boost/shared_ptr.hpp>
#include <string>

#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "XdmfArrayType.hpp"
#include "hdf5.h"

static unsigned short int const compression_factor = 9;

hid_t XdmfType2Hdf5Type(boost::shared_ptr<XdmfArrayType const> xdmf)
{
    std::map<boost::shared_ptr<XdmfArrayType const>, hid_t> xdmf2HdfType = {
        {XdmfArrayType::Float64(), H5T_IEEE_F64LE},
        {XdmfArrayType::Float32(), H5T_IEEE_F32LE},
        {XdmfArrayType::Int32(), H5T_STD_I32LE},
        //{XdmfArrayType::Int64(), H5T_STD_I64LE},
        {XdmfArrayType::UInt32(), H5T_STD_U32LE},
        //{XdmfArrayType::UInt64(), H5T_STD_U64LE},
        {XdmfArrayType::Int8(), H5T_STD_I8LE},
        {XdmfArrayType::UInt8(), H5T_STD_U8LE}};
    try
    {
        return xdmf2HdfType.at(xdmf);
    }
    catch (const std::exception& e)
    {
        OGS_FATAL("No known HDF5 type for XDMF type : {:s} .", xdmf->getName());
    }
}
using namespace std::string_literals;

static std::string getTimeSection(int const step)
{
    return "t_"s + std::to_string(step);
}

namespace MeshLib::IO
{
int writeHDF5Step(std::filesystem::path const& filepath, int const step,
                  std::string const& attribute_name, void const* attribute_data,
                  std::vector<unsigned long long> const& attribute_dims,
                  boost::shared_ptr<const XdmfArrayType> data_type,
                  bool const has_compression_lib)
{
    // \TODO (tm) Errhandling, not implemented as we will change from CAPI to
    // C++API negative value is failure
    hid_t file = H5Fopen(filepath.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    std::string const time_section = getTimeSection(step);

    // Open or create Group
    hid_t group_id = [](auto file, auto time_section) {
        if (H5Lexists(file, time_section.c_str(), H5P_DEFAULT) > 0)
        {
            // negative value is failure
            return H5Gopen2(file, time_section.c_str(), H5P_DEFAULT);
        }

        // negative value is failure
        return H5Gcreate2(file, time_section.c_str(), H5P_DEFAULT, H5P_DEFAULT,
                          H5P_DEFAULT);
    }(file, time_section);

    // negative value is failure
    hid_t dataset_property = H5Pcreate(H5P_DATASET_CREATE);
    // Alternativle set H5Pset_layout(dcpl, H5D_CONTIGUOUS);
    herr_t status;
    if (has_compression_lib)
    {
        status = H5Pset_deflate(dataset_property, compression_factor);
    }
    status = H5Pset_chunk(dataset_property, attribute_dims.size(),
                          attribute_dims.data());

    hid_t space =
        H5Screate_simple(attribute_dims.size(), attribute_dims.data(), nullptr);
    hid_t dset = H5Dcreate(group_id, attribute_name.c_str(),
                           XdmfType2Hdf5Type(data_type), space, H5P_DEFAULT,
                           dataset_property, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      attribute_data);

    status = H5Pclose(dataset_property);
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Gclose(group_id);
    status = H5Fclose(file);

    return status;
}

std::pair<int, bool> writeHDF5Initial(
    std::vector<double> const& nodes,
    std::vector<unsigned long long> const& geometry_dims,
    std::vector<int> const& topology,
    int const step,
    std::filesystem::path const& filepath)
{
    // Check if gzip compression is available and can be used for both
    // compression and decompression.
    bool has_compression_lib = true;
    if (htri_t avail = H5Zfilter_avail(H5Z_FILTER_DEFLATE); !avail)
    {
        WARN("gzip filter not available.\n");
        has_compression_lib = false;
    }
    unsigned int filter_info;
    H5Zget_filter_info(H5Z_FILTER_DEFLATE, &filter_info);
    if (!(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ||
        !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED))
    {
        WARN("gzip filter not available for encoding and decoding.\n");
        has_compression_lib = false;
    }

    hid_t dataset_property_geometry = H5Pcreate(H5P_DATASET_CREATE);
    herr_t status;
    // status = H5Pset_layout(dcpl, H5D_CONTIGUOUS);
    if (has_compression_lib)
    {
        status = H5Pset_deflate(dataset_property_geometry, compression_factor);
    }
    status = H5Pset_chunk(dataset_property_geometry, geometry_dims.size(),
                          geometry_dims.data());
    if (status != 0)
    {
        ERR("H5Pset_layout failed for geometry");
    }

    hid_t geometry_space = H5Screate_simple(2, geometry_dims.data(), nullptr);

    hsize_t topology_dims[1] = {topology.size()};
    hid_t topology_space = H5Screate_simple(1, topology_dims, nullptr);

    {
        hid_t file = H5Fcreate(filepath.string().c_str(), H5F_ACC_TRUNC,
                               H5P_DEFAULT, H5P_DEFAULT);
        std::string const time_section = getTimeSection(step);
        hid_t group_id = H5Gcreate2(file, time_section.c_str(), H5P_DEFAULT,
                                    H5P_DEFAULT, H5P_DEFAULT);

        hid_t dataset_geometry =
            H5Dcreate(group_id, "geometry", H5T_IEEE_F64BE, geometry_space,
                      H5P_DEFAULT, dataset_property_geometry, H5P_DEFAULT);

        status = H5Dwrite(dataset_geometry, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, nodes.data());

        hid_t dataset_property_topology = H5Pcreate(H5P_DATASET_CREATE);

        // status = H5Pset_layout(dcpl, H5D_CONTIGUOUS);
        status = H5Pset_deflate(dataset_property_topology, compression_factor);
        status = H5Pset_chunk(dataset_property_topology, 1, topology_dims);

        hid_t dataset_topology =
            H5Dcreate(group_id, "topology", H5T_STD_I32LE, topology_space,
                      H5P_DEFAULT, dataset_property_topology, H5P_DEFAULT);

        status = H5Dwrite(dataset_topology, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, topology.data());
        status = H5Pclose(dataset_property_geometry);
        status = H5Pclose(dataset_property_topology);
        status = H5Dclose(dataset_geometry);
        status = H5Dclose(dataset_topology);
        status = H5Sclose(geometry_space);
        status = H5Gclose(group_id);
        status = H5Fclose(file);
    }

    return {0, has_compression_lib};
}
}  // namespace MeshLib::IO
