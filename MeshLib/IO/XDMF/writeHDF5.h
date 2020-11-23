/**
 * \file
 * \author Tobias Meisel
 * \date   2020-11-13
 * \brief  Write vectorized data into hdf5 file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <filesystem>
#include <vector>

namespace MeshLib::IO
{
/**
 * \brief Writes the initial data (constant data) to the h5 file.
 * The file will be overwritten if it already exists
 * @param nodes flattened vector of all nodes (points)
 * (containing point coordinates in form  px0 py0 pz0 px1 py1 pz1)
 * @param geometry_dims vector of values at the supporting points
 * @param topology flattened vector of all cells
 * (containg [cell1, cell2] with flattened values for each cell
 * [celltype cellpoint_1..cellpoint_n]
 * @param step number of the timestep
 * @param filepath absolute or relativ filepath to the hdf5 file
 * @return pair with pair.first: returns negative value on error
 * pair.second return true if compression is found on system otherwise
 * false
 */
std::pair<int, bool> writeHDF5Initial(
    std::vector<double> const& nodes,
    std::vector<unsigned long long> const& geometry_dims,
    std::vector<int> const& topology,
    int step,
    std::filesystem::path const& filepath);
/**
 * \brief Appends the data for the time step (constant data) to the h5 file.
 * @param filepath absolute or relativ filepath to the hdf5 file
 * @param step number of the timestep
 * @param attribute_name name of the attribute (property)
 * @param attribute_data pointer to first element of the actual attribute data
 * @param attribute_dims vector of dimensions of the attribute data
 * @param attribue_data_type XdmfArrayType of underlying (e.g. float32)
 * @return returns negative value on error
 */
int writeHDF5Step(std::filesystem::path const& filepath, int step,
                  std::string const& attribute_name, void const* attribute_data,
                  std::vector<unsigned long long> const& attribute_dims,
                  boost::shared_ptr<const XdmfArrayType> attribue_data_type,
                  bool has_compression_lib);
}  // namespace MeshLib::IO
