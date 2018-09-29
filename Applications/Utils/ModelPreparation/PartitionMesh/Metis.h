/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <vector>

namespace MeshLib
{
class Element;
}

namespace ApplicationUtils
{
/// Write elements as METIS graph file
/// \param elements The mesh elements.
/// \param file_name File name with an extension of mesh.
void writeMETIS(std::vector<MeshLib::Element*> const& elements,
                const std::string& file_name);

/// Read metis data
/// \param file_name_base The prefix of the filename.
/// \param number_of_partitions The number is used to compose the full filename
///                             and forms the postfix.
/// \param number_of_nodes Expected/required number of nodes to be read.
std::vector<std::size_t> readMetisData(const std::string& file_name_base,
                                       long number_of_partitions,
                                       std::size_t number_of_nodes);

/// Removes the F.mesh.npart.P and F.mesh.epart.P files, where F is file name
/// base and P is the number of partitions.
void removeMetisPartitioningFiles(std::string const& file_name_base,
                                  long number_of_partitions);

}  // namespace ApplicationUtils
