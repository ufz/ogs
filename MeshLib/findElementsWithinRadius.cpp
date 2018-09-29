/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "findElementsWithinRadius.h"


#include <vector>
#include <unordered_set>

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

    // Collect start element node coordinates into local contigiuos memory.
    std::vector<MathLib::Point3d> start_element_nodes;
    {
        auto const start_element_n_nodes = start_element.getNumberOfNodes();
        start_element_nodes.reserve(start_element_n_nodes);
        for (unsigned n = 0; n < start_element_n_nodes; ++n)
            start_element_nodes.push_back(*start_element.getNode(n));
    }

    // Returns true if the test node is inside the radius of any of the
    // element's nodes.
    auto node_inside_radius = [&start_element_nodes,
                               radius_squared](Node const* test_node) {
        for (auto const& n : start_element_nodes)
        {
            if (MathLib::sqrDist(*test_node, n) <= radius_squared)
                return true;
        }
        return false;
    };

    // Returns true if any of the test element's nodes is inside the start
    // element's radius.
    auto element_in_radius = [&node_inside_radius](Element const& element) {
        auto const n_nodes = element.getNumberOfNodes();
        for (unsigned i = 0; i < n_nodes; ++i)
        {
            if (node_inside_radius(element.getNode(i)))
                return true;
        }
        return false;
    };

    std::set<std::size_t> found_elements;
    std::vector<Element const*> neighbors_to_visit;
    std::unordered_set<std::size_t> visited_elements;

    auto mark_visited_and_add_neighbors =
        [&found_elements, &neighbors_to_visit, &visited_elements](
            Element const& element) {
            found_elements.insert(element.getID());
            visited_elements.insert(element.getID());
            auto const n_neighbors = element.getNumberOfNeighbors();
            for (unsigned n = 0; n < n_neighbors; ++n)
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
        auto const& current_element = *neighbors_to_visit.back();
        neighbors_to_visit.pop_back();

        // If any node is inside the radius, all neighbors are visited.
        if (element_in_radius(current_element))
            mark_visited_and_add_neighbors(current_element);
    }

    return {std::begin(found_elements), std::end(found_elements)};
}
}
