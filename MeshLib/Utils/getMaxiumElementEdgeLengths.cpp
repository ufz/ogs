/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 2, 2022, 3:05 PM
 */

#include "getMaxiumElementEdgeLengths.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
std::vector<double> getMaxiumElementEdgeLengths(
    std::vector<Element*> const& elements)
{
    std::vector<double> element_edge_lengths(elements.size());
    for (auto const element : elements)
    {
        assert(element->getGeomType() != MeshElemType::POINT);

        auto const& [min_edge_length, max_edge_length] =
            MeshLib::computeSqrEdgeLengthRange(*element);
        (void)min_edge_length;

        element_edge_lengths[element->getID()] = std::sqrt(max_edge_length);
    }

    return element_edge_lengths;
}
}  // namespace MeshLib
