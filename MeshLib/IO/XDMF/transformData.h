/**
 * \file
 * \author Tobias Meisel
 * \date   2020-11-13
 * \brief  Transforms OGS Mesh into vectorized data.
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <set>
#include <string>

#include "XdmfHdfData.h"

namespace MeshLib
{
class Mesh;
}

namespace MeshLib::IO
{
/**
 * \brief Create meta data for attributes used for hdf5 and xdmf
 * \param mesh OGS mesh can be mesh or partitionedMesh
 * \param n_files specifies the number of files. If greater than 1 it groups the
 * data of each process to n_files
 * @param chunk_size_bytes Data will be split into chunks. The parameter
 * specifies the size (in bytes) of the largest chunk.
 * @return vector of meta data
 */
std::vector<XdmfHdfData> transformAttributes(MeshLib::Mesh const& mesh,
                                             unsigned int n_files,
                                             unsigned int chunk_size_bytes);
/**
 * \brief Create meta data for geometry used for hdf5 and xdmf
 * \param mesh OGS mesh can be mesh or partitionedMesh
 * \param data_ptr Memory location of geometry values.
 * \param n_files specifies the number of files. If greater than 1 it groups the
 * data of each process to n_files
 * \param chunk_size_bytes Data will be split into chunks. The parameter
 * specifies the size (in bytes) of the largest chunk.
 * \return Geometry meta data
 */
XdmfHdfData transformGeometry(MeshLib::Mesh const& mesh, double const* data_ptr,
                              unsigned int n_files,
                              unsigned int chunk_size_bytes);
/**
 * \brief  Create meta data for topology used for HDF5 and XDMF
 * \param values actual topology values to get size and memory location
 * \param n_files specifies the number of files. If greater than 1 it groups the
 * data of each process to n_files
 * \param chunk_size_bytes Data will be split into chunks. The parameter
 * specifies the size (in bytes) of the largest chunk.
 * \return Topology meta data
 */
XdmfHdfData transformTopology(std::vector<int> const& values,
                              unsigned int n_files,
                              unsigned int chunk_size_bytes);
/**
 * \brief Copies all node points into a new vector. Contiguous data used for
 * writing. Conform with XDMF standard in
 * https://xdmf.org/index.php/XDMF_Model_and_Format
 * \param mesh OGS mesh can be mesh or partitionedMesh
 * \return vector containing a copy of the data
 */
std::vector<double> transformToXDMFGeometry(MeshLib::Mesh const& mesh);
/**
 * \brief Copies all cells into a new vector. Contiguous data used for writing.
 * The topology is specific to xdmf because it contains the xdmf cell types!!
 * See section topology in https://xdmf.org/index.php/XDMF_Model_and_Format
 * \param mesh OGS mesh can be mesh or partitionedMesh
 * \param offset Local offset to transform local to global cell ID. Offset must
 * be zero in serial and must be defined for each process in parallel execution.
 * \return vector containing a copy of the data
 */
std::vector<int> transformToXDMFTopology(MeshLib::Mesh const& mesh,
                                         std::size_t const offset);
}  // namespace MeshLib::IO
