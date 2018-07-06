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

#include <iostream>

#include "BaseLib/Error.h"
#include "MeshLib/Elements/Element.h"

namespace ApplicationUtils
{
void writeMETIS(std::vector<MeshLib::Element*> const& elements,
                const std::string& file_name)
{
    std::ofstream os(file_name, std::ios::trunc);
    if (!os.is_open())
    {
        OGS_FATAL("Error: cannot open file %s.", file_name.data());
    }

    if (!os.good())
    {
        OGS_FATAL("Error: Cannot write in file %s.", file_name.data());
    }

    os << elements.size() << " \n";
    for (const auto* elem : elements)
    {
        os << elem->getNodeIndex(0) + 1;
        for (unsigned j = 1; j < elem->getNumberOfNodes(); j++)
        {
            os << " " << elem->getNodeIndex(j) + 1;
        }
        os << "\n";
    }
}

std::vector<std::size_t> readMetisData(const std::string& file_name_base,
                                       long const number_of_partitions,
                                       std::size_t const number_of_nodes)
{
    const std::string npartitions_str = std::to_string(number_of_partitions);

    // Read partitioned mesh data from METIS
    const std::string fname_parts =
        file_name_base + ".mesh.npart." + npartitions_str;

    std::ifstream npart_in(fname_parts);
    if (!npart_in.is_open())
    {
        OGS_FATAL(
            "Error: cannot open file %s. It may not exist!\n"
            "Run mpmetis beforehand or use option -m",
            fname_parts.data());
    }

    std::vector<std::size_t> partition_ids(number_of_nodes);


    std::size_t counter = 0;
    while (!npart_in.eof())
    {
        npart_in >> partition_ids[counter++] >> std::ws;
        if (counter == number_of_nodes)
        {
            break;
        }
    }

    if (npart_in.bad())
    {
        OGS_FATAL("Error while reading file %s.", fname_parts.data());
    }

    if (counter != number_of_nodes)
    {
        OGS_FATAL("Error: data in %s are less than expected.",
                  fname_parts.data());
    }

    return partition_ids;
}

void removeMetisPartitioningFiles(std::string const& file_name_base,
                                  long const number_of_partitions)
{
    const std::string npartitions_str = std::to_string(number_of_partitions);

    std::remove((file_name_base + ".mesh.npart." + npartitions_str).c_str());
    std::remove((file_name_base + ".mesh.epart." + npartitions_str).c_str());
}
}  // namespace ApplicationUtils
