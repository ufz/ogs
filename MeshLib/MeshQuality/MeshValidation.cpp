/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of the MeshValidation class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshValidation.h"

#include <algorithm>
#include <numeric>
#include <stack>

#include "BaseLib/Logging.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/MeshRevision.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/MeshSurfaceExtraction.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
/** Finds all surface elements that can be reached from element. All elements
 * that are found in this way are marked in the global sfc_idx vector using the
 * current_index.
 * @param element The mesh element from which the search is started
 * @param sfc_idx The global index vector notifying to which surface elements
 * belong
 * @param current_index The index that all elements reachable from element will
 * be assigned in sfc_idx
 */
static void trackSurface(MeshLib::Element const* const element,
                         std::vector<unsigned>& sfc_idx,
                         unsigned const current_index)
{
    std::stack<MeshLib::Element const*> elem_stack;
    elem_stack.push(element);
    while (!elem_stack.empty())
    {
        MeshLib::Element const* const elem = elem_stack.top();
        elem_stack.pop();
        sfc_idx[elem->getID()] = current_index;
        std::size_t const n_neighbors(elem->getNumberOfNeighbors());
        for (std::size_t i = 0; i < n_neighbors; ++i)
        {
            MeshLib::Element const* neighbor(elem->getNeighbor(i));
            if (neighbor != nullptr && sfc_idx[neighbor->getID()] ==
                                           std::numeric_limits<unsigned>::max())
            {
                elem_stack.push(neighbor);
            }
        }
    }
}

MeshValidation::MeshValidation(MeshLib::Mesh& mesh)
{
    INFO("Mesh Quality Control:");
    INFO("Looking for unused nodes...");
    NodeSearch ns(mesh);
    ns.searchUnused();
    if (!ns.getSearchedNodeIDs().empty())
    {
        INFO("{:d} unused mesh nodes found.", ns.getSearchedNodeIDs().size());
    }
    MeshRevision rev(mesh);
    INFO("Found {:d} potentially collapsible nodes.",
         rev.getNumberOfCollapsableNodes());

    const std::vector<ElementErrorCode> codes(
        MeshLib::MeshValidation::testElementGeometry(mesh));
    std::array<std::string,
               static_cast<std::size_t>(ElementErrorFlag::MaxValue)>
        output_str(MeshLib::MeshValidation::ElementErrorCodeOutput(codes));
    for (auto& i : output_str)
    {
        INFO("{:s}", i);
    }
}

std::vector<ElementErrorCode> MeshValidation::testElementGeometry(
    const MeshLib::Mesh& mesh, double min_volume)
{
    INFO("Testing mesh element geometry:");
    const auto nErrorCodes(
        static_cast<std::size_t>(ElementErrorFlag::MaxValue));
    unsigned error_count[nErrorCodes];
    std::fill_n(error_count, 4, 0);
    const std::size_t nElements(mesh.getNumberOfElements());
    const std::vector<MeshLib::Element*>& elements(mesh.getElements());
    std::vector<ElementErrorCode> error_code_vector;
    error_code_vector.reserve(nElements);

    for (std::size_t i = 0; i < nElements; ++i)
    {
        const ElementErrorCode e = elements[i]->validate();
        error_code_vector.push_back(e);
        if (e.none())
        {
            continue;
        }

        // increment error statistics
        const std::bitset<static_cast<std::size_t>(ElementErrorFlag::MaxValue)>
            flags(static_cast<std::bitset<static_cast<std::size_t>(
                      ElementErrorFlag::MaxValue)>>(e));
        for (unsigned j = 0; j < nErrorCodes; ++j)
        {
            error_count[j] += flags[j];
        }
    }

    // if a larger volume threshold is given, evaluate elements again to add
    // them even if they are formally okay
    if (min_volume > std::numeric_limits<double>::epsilon())
    {
        for (std::size_t i = 0; i < nElements; ++i)
        {
            if (elements[i]->getContent() < min_volume)
            {
                error_code_vector[i].set(ElementErrorFlag::ZeroVolume);
            }
        }
    }

    // output
    const auto error_sum(static_cast<unsigned>(
        std::accumulate(error_count, error_count + nErrorCodes, 0.0)));
    if (error_sum != 0)
    {
        ElementErrorFlag flags[nErrorCodes] = {
            ElementErrorFlag::ZeroVolume, ElementErrorFlag::NonCoplanar,
            ElementErrorFlag::NonConvex, ElementErrorFlag::NodeOrder};
        for (std::size_t i = 0; i < nErrorCodes; ++i)
        {
            if (error_count[i])
            {
                INFO("{:d} elements found with {:s}.",
                     error_count[i],
                     ElementErrorCode::toString(flags[i]));
            }
        }
    }
    else
    {
        INFO("No errors found.");
    }
    return error_code_vector;
}

std::array<std::string, static_cast<std::size_t>(ElementErrorFlag::MaxValue)>
MeshValidation::ElementErrorCodeOutput(
    const std::vector<ElementErrorCode>& error_codes)
{
    const auto nErrorFlags(
        static_cast<std::size_t>(ElementErrorFlag::MaxValue));
    ElementErrorFlag flags[nErrorFlags] = {
        ElementErrorFlag::ZeroVolume, ElementErrorFlag::NonCoplanar,
        ElementErrorFlag::NonConvex, ElementErrorFlag::NodeOrder};
    const std::size_t nElements(error_codes.size());
    std::array<std::string,
               static_cast<std::size_t>(ElementErrorFlag::MaxValue)>
        output;

    for (std::size_t i = 0; i < nErrorFlags; ++i)
    {
        unsigned count(0);
        std::string elementIdStr;

        for (std::size_t j = 0; j < nElements; ++j)
        {
            if (error_codes[j][flags[i]])
            {
                elementIdStr += (std::to_string(j) + ", ");
                count++;
            }
        }
        const std::string nErrorsStr = (count) ? std::to_string(count) : "No";
        output[i] = (nErrorsStr + " elements found with " +
                     ElementErrorCode::toString(flags[i]) + ".\n");

        if (count)
        {
            output[i] += ("ElementIDs: " + elementIdStr + "\n");
        }
    }
    return output;
}

unsigned MeshValidation::detectHoles(MeshLib::Mesh const& mesh)
{
    if (mesh.getDimension() == 1)
    {
        return 0;
    }

    auto boundary_mesh = MeshLib::BoundaryExtraction::getBoundaryElementsAsMesh(
        mesh, "bulk_node_ids", "bulk_element_ids", "bulk_face_ids");
    std::vector<MeshLib::Element*> const& elements(
        boundary_mesh->getElements());

    std::vector<unsigned> sfc_idx(elements.size(),
                                  std::numeric_limits<unsigned>::max());
    unsigned current_surface_id(0);
    auto it = sfc_idx.cbegin();

    while (it != sfc_idx.cend())
    {
        std::size_t const idx =
            static_cast<std::size_t>(std::distance(sfc_idx.cbegin(), it));
        trackSurface(elements[idx], sfc_idx, current_surface_id++);
        it = std::find(sfc_idx.cbegin(),
                       sfc_idx.cend(),
                       std::numeric_limits<unsigned>::max());
    }

    // Subtract "1" from the number of surfaces found to get the number of
    // holes.
    return (--current_surface_id);
}
}  // end namespace MeshLib
