/**
 * \file
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "findElementsWithinRadius.h"


#include <deque>
#include <unordered_set>

#include "BaseLib/uniqueInsert.h"
#include "MathLib/MathTools.h"

#include "Elements/Element.h"
#include "Node.h"

namespace MeshLib {
std::vector<std::size_t> findElementsWithinRadius(Element const& start_element,
                                                  double const radius_squared)
{
    // Special case for 0 radius. All other radii will include at least one
    // neighbor of the start element.
    if (radius_squared == 0.)
        return {start_element.getID()};

    // Returns true if the test node is inside the radius of any of the
    // element's nodes.
    auto node_inside_radius = [&start_element,
                               radius_squared](Node const* test_node) {
        for (unsigned n = 0; n < start_element.getNumberOfNodes(); ++n)
        {
            if (MathLib::sqrDist(*test_node, *start_element.getNode(n)) <=
                radius_squared)
                return true;
        }
        return false;
    };

    std::set<std::size_t> found_elements;
    std::deque<Element const*> neighbors_to_visit;
    std::unordered_set<std::size_t> visited_elements;

    auto mark_visited_and_add_neighbors =
        [&found_elements, &neighbors_to_visit,
         &visited_elements](Element const& element) {
            found_elements.insert(element.getID());
            visited_elements.insert(element.getID());
            for (unsigned n = 0; n < element.getNumberOfNeighbors(); ++n)
            {
                auto neighbor = element.getNeighbor(n);
                if (neighbor == nullptr)
                    continue;
                auto x = visited_elements.find(neighbor->getID());
                if (x != visited_elements.end())
                    continue;

                neighbors_to_visit.push_back(neighbor);
                visited_elements.insert(neighbor->getID());
            }
        };

    // The result always contains the starting element.
    mark_visited_and_add_neighbors(start_element);

    while (!neighbors_to_visit.empty())
    {
        auto const& current_element = *neighbors_to_visit.front();
        neighbors_to_visit.pop_front();

        // if any node is inside the radius, all neighbors are visited.
        for (unsigned i = 0; i < current_element.getNumberOfNodes(); ++i)
        {
            if (!node_inside_radius(current_element.getNode(i)))
                continue;

            mark_visited_and_add_neighbors(current_element);
        }
    }

    return {std::begin(found_elements), std::end(found_elements)};
}
}
