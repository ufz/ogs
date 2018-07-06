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

}  // namespace ApplicationUtils
