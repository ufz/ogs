/**
 * \file
 * \author Karsten Rink
 * \date   2014-02-14
 * \brief  Implementation of the MeshRevision class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshRevision.h"

#include <numeric>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>

#include "BaseLib/Algorithm.h"
#include "BaseLib/Logging.h"
#include "GeoLib/Grid.h"
#include "MathLib/GeometricBasics.h"
#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Properties.h"
#include "MeshLib/Utils/DuplicateMeshComponents.h"

namespace MeshToolsLib
{
/// Lookup-table for returning the third node of bottom or top triangle given
/// the other two.
unsigned lutPrismThirdNode(unsigned const id1, unsigned const id2)
{
    if ((id1 == 0 && id2 == 1) || (id1 == 1 && id2 == 0))
    {
        return 2;
    }
    if ((id1 == 1 && id2 == 2) || (id1 == 2 && id2 == 1))
    {
        return 0;
    }
    if ((id1 == 0 && id2 == 2) || (id1 == 2 && id2 == 0))
    {
        return 1;
    }
    if ((id1 == 3 && id2 == 4) || (id1 == 4 && id2 == 3))
    {
        return 5;
    }
    if ((id1 == 4 && id2 == 5) || (id1 == 5 && id2 == 4))
    {
        return 3;
    }
    if ((id1 == 3 && id2 == 5) || (id1 == 5 && id2 == 3))
    {
        return 4;
    }
    return std::numeric_limits<unsigned>::max();
}
}  // namespace MeshToolsLib

namespace
{
template <typename ElementType>
std::unique_ptr<MeshLib::Element> createElement(
    std::span<MeshLib::Node* const> const element_nodes,
    std::vector<MeshLib::Node*> const& nodes,
    std::array<std::size_t, ElementType::n_all_nodes> const local_ids)
{
    using namespace MeshLib::views;
    auto lookup_in = [](auto const& values)
    {
        return ranges::views::transform([&values](std::size_t const n)
                                        { return values[n]; });
    };

    std::array<MeshLib::Node*, ElementType::n_all_nodes> new_nodes;
    ranges::copy(local_ids | lookup_in(element_nodes) | ids | lookup_in(nodes),
                 begin(new_nodes));

    return std::make_unique<ElementType>(new_nodes);
}

/// Subdivides a nonplanar quad into two triangles.
unsigned subdivideQuad(MeshLib::Element const* const quad,
                       std::vector<MeshLib::Node*> const& nodes,
                       std::vector<MeshLib::Element*>& new_elements)
{
    std::array<std::size_t, 3> const tri1_node_ids{0, 1, 2};
    new_elements.push_back(
        createElement<MeshLib::Tri>(quad->nodes(), nodes, tri1_node_ids)
            .release());

    std::array<std::size_t, 3> const tri2_node_ids{0, 2, 3};
    new_elements.push_back(
        createElement<MeshLib::Tri>(quad->nodes(), nodes, tri2_node_ids)
            .release());

    return 2;
}

/// Subdivides a prism with nonplanar quad faces into two tets.
unsigned subdividePrism(MeshLib::Element const* const prism,
                        std::vector<MeshLib::Node*> const& nodes,
                        std::vector<MeshLib::Element*>& new_elements)
{
    auto addTetrahedron =
        [&prism, &nodes, &new_elements](std::array<std::size_t, 4> const ids)
    {
        new_elements.push_back(
            createElement<MeshLib::Tet>(prism->nodes(), nodes, ids).release());
    };

    addTetrahedron({0, 1, 2, 3});
    addTetrahedron({3, 2, 4, 5});
    addTetrahedron({2, 1, 3, 4});

    return 3;
}

/// Subdivides a Hex with nonplanar faces into tets.
unsigned subdivideHex(MeshLib::Element const* const hex,
                      std::vector<MeshLib::Node*> const& nodes,
                      std::vector<MeshLib::Element*>& new_elements)
{
    auto prism1 =
        createElement<MeshLib::Prism>(hex->nodes(), nodes, {0, 2, 1, 4, 6, 5});
    subdividePrism(prism1.get(), nodes, new_elements);

    auto prism2 =
        createElement<MeshLib::Prism>(hex->nodes(), nodes, {4, 6, 7, 0, 2, 3});
    subdividePrism(prism2.get(), nodes, new_elements);

    return 6;
}

/// Subdivides a pyramid with a nonplanar base into two tets.
unsigned subdividePyramid(MeshLib::Element const* const pyramid,
                          std::vector<MeshLib::Node*> const& nodes,
                          std::vector<MeshLib::Element*>& new_elements)
{
    auto addTetrahedron =
        [&pyramid, &nodes, &new_elements](std::array<std::size_t, 4> const ids)
    {
        new_elements.push_back(
            createElement<MeshLib::Tet>(pyramid->nodes(), nodes, ids)
                .release());
    };

    addTetrahedron({0, 1, 2, 4});
    addTetrahedron({0, 2, 3, 4});

    return 2;
}

/// Creates a line element from the first two unique nodes found in the element
/// (element *should* have exactly two unique nodes!)
MeshLib::Element* constructLine(MeshLib::Element const* const element,
                                const std::vector<MeshLib::Node*>& nodes)
{
    std::array<std::size_t, 2> line_node_ids = {0, 0};
    for (unsigned i = 1; i < element->getNumberOfBaseNodes(); ++i)
    {
        if (element->getNode(i)->getID() != element->getNode(0)->getID())
        {
            line_node_ids[1] = i;
            break;
        }
    }
    assert(line_node_ids[1] != 0);
    return createElement<MeshLib::Line>(element->nodes(), nodes, line_node_ids)
        .release();
}

/// Creates a triangle element from the first three unique nodes found in the
/// element (element *should* have exactly three unique nodes!)
MeshLib::Element* constructTri(MeshLib::Element const* const element,
                               const std::vector<MeshLib::Node*>& nodes)
{
    // TODO?
    // In theory three unique nodes could also be reduced to two lines e.g. with
    // a quad where two diametral nodes collapse. This case is currently not
    // implemented!
    std::array<MeshLib::Node*, 3> tri_nodes;
    tri_nodes[0] = nodes[element->getNode(0)->getID()];
    tri_nodes[2] = nullptr;
    for (unsigned i = 1; i < element->getNumberOfBaseNodes(); ++i)
    {
        if (element->getNode(i)->getID() != tri_nodes[0]->getID())
        {
            tri_nodes[1] = nodes[element->getNode(i)->getID()];
            for (unsigned j = i + 1; j < element->getNumberOfBaseNodes(); ++j)
            {
                if (element->getNode(j)->getID() != tri_nodes[1]->getID())
                {
                    tri_nodes[2] = nodes[element->getNode(j)->getID()];
                    break;
                }
            }
            if (tri_nodes[2])
            {
                break;
            }
        }
    }
    assert(tri_nodes[2] != nullptr);
    return new MeshLib::Tri(tri_nodes);
}

/// Creates a quad or a tet, depending if the four nodes being coplanar or not
/// (element *should* have exactly four unique nodes!)
MeshLib::Element* constructFourNodeElement(
    MeshLib::Element const* const element,
    std::vector<MeshLib::Node*> const& nodes,
    unsigned const min_elem_dim = 1)
{
    std::array<MeshLib::Node*, 4> new_nodes;
    unsigned count(0);
    new_nodes[count++] = nodes[element->getNode(0)->getID()];
    for (unsigned i = 1; i < element->getNumberOfBaseNodes(); ++i)
    {
        if (count > 3)
        {
            break;
        }
        bool unique_node(true);
        for (unsigned j = 0; j < i; ++j)
        {
            if (element->getNode(i)->getID() == element->getNode(j)->getID())
            {
                unique_node = false;
                break;
            }
        }
        if (unique_node)
        {
            new_nodes[count++] = nodes[element->getNode(i)->getID()];
        };
    }

    // test if quad or tet
    const bool isQuad(MathLib::isCoplanar(*new_nodes[0], *new_nodes[1],
                                          *new_nodes[2], *new_nodes[3]));
    if (isQuad && min_elem_dim < 3)
    {
        MeshLib::Element* elem(new MeshLib::Quad(new_nodes));
        for (unsigned i = 1; i < 3; ++i)
        {
            if (elem->validate().none())
            {
                return elem;
            }

            // change node order if not convex
            MeshLib::Node* tmp = new_nodes[i + 1];
            new_nodes[i + 1] = new_nodes[i];
            new_nodes[i] = tmp;
        }
        return elem;
    }
    if (!isQuad)
    {
        return new MeshLib::Tet(new_nodes);
    }
    // is quad but min elem dim == 3

    return nullptr;
}

/// Reduces a pyramid element by removing collapsed nodes and constructing a new
/// elements from the remaining nodes.
void reducePyramid(MeshLib::Element const* const org_elem,
                   unsigned const n_unique_nodes,
                   std::vector<MeshLib::Node*> const& nodes,
                   std::vector<MeshLib::Element*>& new_elements,
                   unsigned const min_elem_dim)
{
    if (n_unique_nodes == 4)
    {
        MeshLib::Element* elem(
            constructFourNodeElement(org_elem, nodes, min_elem_dim));
        if (elem)
        {
            new_elements.push_back(elem);
        }
    }
    else if (n_unique_nodes == 3 && min_elem_dim < 3)
    {
        new_elements.push_back(constructTri(org_elem, nodes));
    }
    else if (n_unique_nodes == 2 && min_elem_dim == 1)
    {
        new_elements.push_back(constructLine(org_elem, nodes));
    }
}

/// Reduces a prism element by removing collapsed nodes and constructing one or
/// two new elements from the remaining nodes.
/// @return The number of newly created elements
unsigned reducePrism(MeshLib::Element const* const org_elem,
                     unsigned const n_unique_nodes,
                     std::vector<MeshLib::Node*> const& nodes,
                     std::vector<MeshLib::Element*>& new_elements,
                     unsigned const min_elem_dim)
{
    auto addTetrahedron =
        [&org_elem, &nodes, &new_elements](std::array<std::size_t, 4> const ids)
    {
        new_elements.push_back(
            createElement<MeshLib::Tet>(org_elem->nodes(), nodes, ids)
                .release());
    };

    // TODO?
    // In theory a node from the bottom triangle and a node from the top
    // triangle that are not connected by an edge could collapse, resulting in a
    // combination of tri and quad elements. This case is currently not tested.

    // if one of the non-triangle edges collapsed, elem can be reduced to a
    // pyramid, otherwise it will be two tets
    if (n_unique_nodes == 5)
    {
        for (unsigned i = 0; i < 5; ++i)
        {
            for (unsigned j = i + 1; j < 6; ++j)
            {
                if (i != j && org_elem->getNode(i)->getID() ==
                                  org_elem->getNode(j)->getID())
                {
                    // non triangle edge collapsed
                    if (i % 3 == j % 3)
                    {
                        addTetrahedron(
                            {(i + 1) % 3, (i + 2) % 3, i, (i + 1) % 3 + 3});
                        addTetrahedron(
                            {(i + 1) % 3 + 3, (i + 2) % 3, i, (i + 2) % 3 + 3});
                        return 2;
                    }

                    // triangle edge collapsed
                    const unsigned i_offset = (i > 2) ? i - 3 : i + 3;
                    const unsigned j_offset = (i > 2) ? j - 3 : j + 3;
                    const unsigned k = MeshToolsLib::lutPrismThirdNode(i, j);
                    if (k == std::numeric_limits<unsigned>::max())
                    {
                        ERR("Unexpected error during prism reduction.");
                        return 0;
                    }
                    const unsigned k_offset = (i > 2) ? k - 3 : k + 3;

                    addTetrahedron({i_offset, j_offset, k_offset, i});

                    const unsigned l =
                        (MathLib::isCoplanar(*org_elem->getNode(i_offset),
                                             *org_elem->getNode(k_offset),
                                             *org_elem->getNode(i),
                                             *org_elem->getNode(k)))
                            ? j
                            : i;
                    const unsigned l_offset = (i > 2) ? l - 3 : l + 3;
                    addTetrahedron({l_offset, k_offset, i, k});
                    return 2;
                }
            }
        }
    }
    else if (n_unique_nodes == 4)
    {
        MeshLib::Element* const elem(
            constructFourNodeElement(org_elem, nodes, min_elem_dim));
        if (elem)
        {
            new_elements.push_back(elem);
        }
    }
    else if (n_unique_nodes == 3 && min_elem_dim < 3)
    {
        new_elements.push_back(constructTri(org_elem, nodes));
    }
    else if (n_unique_nodes == 2 && min_elem_dim == 1)
    {
        new_elements.push_back(constructLine(org_elem, nodes));
    }
    return 1;
}

/// Lookup-table for returning four nodes connected to the two nodes (id1, id2)
/// forming an edge in a Hex.
std::array<unsigned, 4> lutHexCuttingQuadNodes(unsigned id1, unsigned id2)
{
    if (id1 == 0 && id2 == 1)
    {
        return {3, 2, 5, 4};
    }
    if (id1 == 1 && id2 == 2)
    {
        return {0, 3, 6, 5};
    }
    if (id1 == 2 && id2 == 3)
    {
        return {1, 0, 7, 6};
    }
    if (id1 == 3 && id2 == 0)
    {
        return {2, 1, 4, 7};
    }
    if (id1 == 4 && id2 == 5)
    {
        return {0, 1, 6, 7};
    }
    if (id1 == 5 && id2 == 6)
    {
        return {1, 2, 7, 4};
    }
    if (id1 == 6 && id2 == 7)
    {
        return {2, 3, 4, 5};
    }
    if (id1 == 7 && id2 == 4)
    {
        return {3, 0, 5, 6};
    }
    if (id1 == 0 && id2 == 4)
    {
        return {3, 7, 5, 1};
    }
    if (id1 == 1 && id2 == 5)
    {
        return {0, 4, 6, 2};
    }
    if (id1 == 2 && id2 == 6)
    {
        return {1, 5, 7, 3};
    }
    if (id1 == 3 && id2 == 7)
    {
        return {2, 6, 4, 0};
    }
    if (id1 == 1 && id2 == 0)
    {
        return {2, 3, 4, 5};
    }
    if (id1 == 2 && id2 == 1)
    {
        return {3, 0, 5, 6};
    }
    if (id1 == 3 && id2 == 2)
    {
        return {0, 1, 6, 7};
    }
    if (id1 == 0 && id2 == 3)
    {
        return {1, 2, 7, 4};
    }
    if (id1 == 5 && id2 == 4)
    {
        return {1, 0, 7, 6};
    }
    if (id1 == 6 && id2 == 5)
    {
        return {2, 1, 4, 7};
    }
    if (id1 == 7 && id2 == 6)
    {
        return {3, 2, 5, 4};
    }
    if (id1 == 4 && id2 == 7)
    {
        return {0, 3, 6, 5};
    }
    if (id1 == 4 && id2 == 0)
    {
        return {7, 3, 1, 5};
    }
    if (id1 == 5 && id2 == 1)
    {
        return {4, 0, 2, 6};
    }
    if (id1 == 6 && id2 == 2)
    {
        return {5, 1, 3, 7};
    }
    if (id1 == 7 && id2 == 3)
    {
        return {6, 2, 0, 4};
    }

    OGS_FATAL(
        "lutHexCuttingQuadNodes() for nodes {} and {} does not have a valid "
        "return value.",
        id1, id2);
}

/// Lookup-table for returning the diametral node id of the given node id in a
/// Hex.
unsigned lutHexDiametralNode(unsigned const id)
{
    constexpr std::array<unsigned, 8> hex_diametral_node_ids = {
        {6, 7, 4, 5, 2, 3, 0, 1}};

    return hex_diametral_node_ids[id];
}

/// When a hex is subdivided into two prisms, this returns the nodes of the hex
/// edge that will serve as the back of one of the prisms.
std::pair<unsigned, unsigned> lutHexBackNodes(unsigned const i,
                                              unsigned const j,
                                              unsigned const k,
                                              unsigned const l)
{
    // collapsed edges are *not* connected
    if (lutHexDiametralNode(i) == k)
    {
        return {i, lutHexDiametralNode(l)};
    }
    if (lutHexDiametralNode(i) == l)
    {
        return {i, lutHexDiametralNode(k)};
    }
    if (lutHexDiametralNode(j) == k)
    {
        return {j, lutHexDiametralNode(l)};
    }
    if (lutHexDiametralNode(j) == l)
    {
        return {j, lutHexDiametralNode(k)};
    }

    // collapsed edges *are* connected
    if (i == k)
    {
        return {lutHexDiametralNode(l), j};
    }
    if (i == l)
    {
        return {lutHexDiametralNode(k), j};
    }
    if (j == k)
    {
        return {lutHexDiametralNode(l), i};
    }
    if (j == l)
    {
        return {lutHexDiametralNode(k), i};
    }

    return {std::numeric_limits<unsigned>::max(),
            std::numeric_limits<unsigned>::max()};
}

// In an element with 5 unique nodes, return the node that will be the top of
// the resulting pyramid.
unsigned findPyramidTopNode(MeshLib::Element const& element,
                            std::array<std::size_t, 4> const& base_node_ids)
{
    const std::size_t nNodes(element.getNumberOfBaseNodes());
    for (std::size_t i = 0; i < nNodes; ++i)
    {
        bool top_node = true;
        for (unsigned j = 0; j < 4; ++j)
        {
            if (element.getNode(i)->getID() == base_node_ids[j])
            {
                top_node = false;
            }
        }
        if (top_node)
        {
            return i;
        }
    }
    return std::numeric_limits<unsigned>::max();  // should never be reached if
                                                  // called correctly
}

/// Reduces a hexahedron element by removing collapsed nodes and constructing
/// one or more new elements from the remaining nodes.
/// @return The number of newly created elements
unsigned reduceHex(MeshLib::Element const* const org_elem,
                   unsigned const n_unique_nodes,
                   std::vector<MeshLib::Node*> const& nodes,
                   std::vector<MeshLib::Element*>& new_elements,
                   unsigned const min_elem_dim)
{
    // TODO?
    // if two diametral nodes collapse, all kinds of bizarre (2D-)element
    // combinations could be the result. this case is currently not implemented!

    if (n_unique_nodes == 7)
    {
        // reduce to prism + pyramid
        for (unsigned i = 0; i < 7; ++i)
        {
            for (unsigned j = i + 1; j < 8; ++j)
            {
                if (org_elem->getNode(i)->getID() ==
                    org_elem->getNode(j)->getID())
                {
                    const std::array<unsigned, 4> base_node_ids(
                        lutHexCuttingQuadNodes(i, j));
                    std::array<std::size_t, 5> const pyr_node_ids = {
                        base_node_ids[0], base_node_ids[1], base_node_ids[2],
                        base_node_ids[3], i};
                    new_elements.push_back(
                        createElement<MeshLib::Pyramid>(org_elem->nodes(),
                                                        nodes, pyr_node_ids)
                            .release());

                    if (i < 4 && j >= 4)
                    {
                        std::swap(i, j);
                    }
                    std::array<std::size_t, 6> const prism_node_ids{
                        base_node_ids[0],       base_node_ids[3],
                        lutHexDiametralNode(j), base_node_ids[1],
                        base_node_ids[2],       lutHexDiametralNode(i)};
                    new_elements.push_back(
                        createElement<MeshLib::Prism>(org_elem->nodes(), nodes,
                                                      prism_node_ids)
                            .release());
                    return 2;
                }
            }
        }
    }
    else if (n_unique_nodes == 6)
    {
        // reduce to prism
        for (unsigned i = 0; i < 6; ++i)
        {
            const MeshLib::Element* face(org_elem->getFace(i));
            if (face->getNode(0)->getID() == face->getNode(1)->getID() &&
                face->getNode(2)->getID() == face->getNode(3)->getID())
            {
                std::array<std::size_t, 6> const prism_node_ids{
                    lutHexDiametralNode(
                        getNodeIDinElement(*org_elem, face->getNode(0))),
                    lutHexDiametralNode(
                        getNodeIDinElement(*org_elem, face->getNode(1))),
                    getNodeIDinElement(*org_elem, face->getNode(2)),
                    lutHexDiametralNode(
                        getNodeIDinElement(*org_elem, face->getNode(2))),
                    lutHexDiametralNode(
                        getNodeIDinElement(*org_elem, face->getNode(3))),
                    getNodeIDinElement(*org_elem, face->getNode(0))};

                new_elements.push_back(
                    createElement<MeshLib::Prism>(org_elem->nodes(), nodes,
                                                  prism_node_ids)
                        .release());
                delete face;
                return 1;
            }
            if (face->getNode(0)->getID() == face->getNode(3)->getID() &&
                face->getNode(1)->getID() == face->getNode(2)->getID())
            {
                std::array<std::size_t, 6> const prism_node_ids{
                    lutHexDiametralNode(
                        getNodeIDinElement(*org_elem, face->getNode(0))),
                    lutHexDiametralNode(
                        getNodeIDinElement(*org_elem, face->getNode(3))),
                    getNodeIDinElement(*org_elem, face->getNode(2)),
                    lutHexDiametralNode(
                        getNodeIDinElement(*org_elem, face->getNode(1))),
                    lutHexDiametralNode(
                        getNodeIDinElement(*org_elem, face->getNode(2))),
                    getNodeIDinElement(*org_elem, face->getNode(0))};
                new_elements.push_back(
                    createElement<MeshLib::Prism>(org_elem->nodes(), nodes,
                                                  prism_node_ids)
                        .release());
                delete face;
                return 1;
            }
            delete face;
        }
        // reduce to four tets -> divide into 2 prisms such that each has one
        // collapsed node
        for (unsigned i = 0; i < 7; ++i)
        {
            for (unsigned j = i + 1; j < 8; ++j)
            {
                if (org_elem->getNode(i)->getID() ==
                    org_elem->getNode(j)->getID())
                {
                    for (unsigned k = i; k < 7; ++k)
                    {
                        for (unsigned l = k + 1; l < 8; ++l)
                        {
                            if (!(i == k && j == l) && org_elem->isEdge(i, j) &&
                                org_elem->isEdge(k, l) &&
                                org_elem->getNode(k)->getID() ==
                                    org_elem->getNode(l)->getID())
                            {
                                const std::pair<unsigned, unsigned> back(
                                    lutHexBackNodes(i, j, k, l));
                                if (back.first ==
                                        std::numeric_limits<unsigned>::max() ||
                                    back.second ==
                                        std::numeric_limits<unsigned>::max())
                                {
                                    ERR("Unexpected error during Hex "
                                        "reduction");
                                    return 0;
                                }

                                std::array<unsigned, 4> const cutting_plane(
                                    lutHexCuttingQuadNodes(back.first,
                                                           back.second));
                                std::array<std::size_t, 6> const pris1_node_ids{
                                    back.first,       cutting_plane[0],
                                    cutting_plane[3], back.second,
                                    cutting_plane[1], cutting_plane[2]};
                                auto prism1 = createElement<MeshLib::Prism>(
                                    org_elem->nodes(), nodes, pris1_node_ids);
                                unsigned nNewElements =
                                    reducePrism(prism1.get(), 5, nodes,
                                                new_elements, min_elem_dim);

                                std::array<std::size_t, 6> const pris2_node_ids{
                                    lutHexDiametralNode(back.first),
                                    cutting_plane[0],
                                    cutting_plane[3],
                                    lutHexDiametralNode(back.second),
                                    cutting_plane[1],
                                    cutting_plane[2]};
                                auto prism2 = createElement<MeshLib::Prism>(
                                    org_elem->nodes(), nodes, pris2_node_ids);
                                nNewElements +=
                                    reducePrism(prism2.get(), 5, nodes,
                                                new_elements, min_elem_dim);
                                return nNewElements;
                            }
                        }
                    }
                }
            }
        }
    }
    else if (n_unique_nodes == 5)
    {
        MeshLib::Element* tet1(constructFourNodeElement(org_elem, nodes));
        std::array<std::size_t, 4> const first_four_node_ids = {
            {tet1->getNode(0)->getID(), tet1->getNode(1)->getID(),
             tet1->getNode(2)->getID(), tet1->getNode(3)->getID()}};
        unsigned const fifth_node =
            findPyramidTopNode(*org_elem, first_four_node_ids);

        bool tet_changed(false);
        if (tet1->getGeomType() == MeshLib::MeshElemType::QUAD)
        {
            delete tet1;
            tet_changed = true;
            std::array const tet1_nodes = {
                nodes[first_four_node_ids[0]], nodes[first_four_node_ids[1]],
                nodes[first_four_node_ids[2]],
                nodes[org_elem->getNode(fifth_node)->getID()]};
            new_elements.push_back(new MeshLib::Tet(tet1_nodes));
        }
        else
        {
            new_elements.push_back(tet1);
        }

        std::array const tet2_nodes = {
            (tet_changed) ? nodes[first_four_node_ids[0]]
                          : nodes[first_four_node_ids[1]],
            nodes[first_four_node_ids[2]], nodes[first_four_node_ids[3]],
            nodes[org_elem->getNode(fifth_node)->getID()]};
        new_elements.push_back(new MeshLib::Tet(tet2_nodes));
        return 2;
    }
    else if (n_unique_nodes == 4)
    {
        MeshLib::Element* elem(
            constructFourNodeElement(org_elem, nodes, min_elem_dim));
        if (elem)
        {
            new_elements.push_back(elem);
            return 1;
        }
    }
    else if (n_unique_nodes == 3 && min_elem_dim < 3)
    {
        new_elements.push_back(constructTri(org_elem, nodes));
        return 1;
    }
    else if (min_elem_dim == 1)
    {
        new_elements.push_back(constructLine(org_elem, nodes));
        return 1;
    }
    return 0;
}

/// Subdivides an element if it has a face that is not coplanar
/// @param element the element that will be subdivided
/// @param nodes vector containing the nodes the elements originated by the
/// subdivision are based on
/// @param elements vector of MeshLib::Elements; the elements originated by the
/// subdivision will be inserted into elements
/// @return the number of elements originated by the subdivision
std::size_t subdivideElement(MeshLib::Element const* const element,
                             std::vector<MeshLib::Node*> const& nodes,
                             std::vector<MeshLib::Element*>& elements)
{
    if (element->getGeomType() == MeshLib::MeshElemType::QUAD)
    {
        return subdivideQuad(element, nodes, elements);
    }
    if (element->getGeomType() == MeshLib::MeshElemType::HEXAHEDRON)
    {
        return subdivideHex(element, nodes, elements);
    }
    if (element->getGeomType() == MeshLib::MeshElemType::PYRAMID)
    {
        return subdividePyramid(element, nodes, elements);
    }
    if (element->getGeomType() == MeshLib::MeshElemType::PRISM)
    {
        return subdividePrism(element, nodes, elements);
    }
    return 0;
}

// Revises an element by removing collapsed nodes, using the nodes vector from
// the result mesh.
std::size_t reduceElement(MeshLib::Element const* const element,
                          unsigned const n_unique_nodes,
                          std::vector<MeshLib::Node*> const& nodes,
                          std::vector<MeshLib::Element*>& elements,
                          unsigned const min_elem_dim)
{
    /***************
     * TODO: modify neighbouring elements if one elements has been subdivided
     ***************/
    if (element->getGeomType() == MeshLib::MeshElemType::TRIANGLE &&
        min_elem_dim == 1)
    {
        elements.push_back(constructLine(element, nodes));
        return 1;
    }
    if ((element->getGeomType() == MeshLib::MeshElemType::QUAD) ||
        (element->getGeomType() == MeshLib::MeshElemType::TETRAHEDRON))
    {
        if (n_unique_nodes == 3 && min_elem_dim < 3)
        {
            elements.push_back(constructTri(element, nodes));
        }
        else if (min_elem_dim == 1)
        {
            elements.push_back(constructLine(element, nodes));
        }
        return 1;
    }
    if (element->getGeomType() == MeshLib::MeshElemType::HEXAHEDRON)
    {
        return reduceHex(element, n_unique_nodes, nodes, elements,
                         min_elem_dim);
    }
    if (element->getGeomType() == MeshLib::MeshElemType::PYRAMID)
    {
        reducePyramid(element, n_unique_nodes, nodes, elements, min_elem_dim);
        return 1;
    }
    if (element->getGeomType() == MeshLib::MeshElemType::PRISM)
    {
        return reducePrism(element, n_unique_nodes, nodes, elements,
                           min_elem_dim);
    }

    ERR("Unknown element type.");
    return 0;
}
/// Calculates the number of unique nodes in an element (i.e. uncollapsed
/// nodes).
unsigned getNumberOfUniqueNodes(MeshLib::Element const* const element)
{
    unsigned const nNodes(element->getNumberOfBaseNodes());
    unsigned count(nNodes);

    for (unsigned i = 0; i < nNodes - 1; ++i)
    {
        for (unsigned j = i + 1; j < nNodes; ++j)
        {
            if (element->getNode(i)->getID() == element->getNode(j)->getID())
            {
                count--;
                break;
            }
        }
    }
    return count;
}

template <typename T>
void fillNodeProperty(std::vector<T>& new_prop,
                      std::vector<T> const& old_prop,
                      std::vector<size_t> const& node_ids)
{
    std::size_t const n_nodes = node_ids.size();
    for (std::size_t i = 0; i < n_nodes; ++i)
    {
        if (node_ids[i] != i)
        {
            continue;
        }
        new_prop.push_back(old_prop[i]);
    }
}

template <typename T>
void fillElemProperty(std::vector<T>& new_prop,
                      std::vector<T> const& old_prop,
                      std::vector<size_t> const& elem_ids)
{
    std::transform(elem_ids.cbegin(), elem_ids.cend(),
                   std::back_inserter(new_prop),
                   [&](std::size_t const i) { return old_prop[i]; });
}

/// Copies all scalar arrays according to the restructured Node- and
/// Element-vectors after the mesh revision process (i.e. collapsed nodes, split
/// elements, etc.)
MeshLib::Properties copyProperties(MeshLib::Properties const& props,
                                   std::vector<std::size_t> const& node_ids,
                                   std::vector<std::size_t> const& elem_ids)
{
    auto const prop_names = props.getPropertyVectorNames();
    MeshLib::Properties new_properties;

    for (auto name : prop_names)
    {
        if (props.existsPropertyVector<int>(name, MeshLib::MeshItemType::Node,
                                            1))
        {
            auto const* p = props.getPropertyVector<int>(
                name, MeshLib::MeshItemType::Node, 1);
            auto new_node_vec = new_properties.createNewPropertyVector<int>(
                name, MeshLib::MeshItemType::Node, 1);
            fillNodeProperty(*new_node_vec, *p, node_ids);
            continue;
        }
        if (props.existsPropertyVector<float>(name, MeshLib::MeshItemType::Node,
                                              1))
        {
            auto const* p = props.getPropertyVector<float>(
                name, MeshLib::MeshItemType::Node, 1);
            auto new_node_vec = new_properties.createNewPropertyVector<float>(
                name, MeshLib::MeshItemType::Node, 1);
            fillNodeProperty(*new_node_vec, *p, node_ids);
            continue;
        }
        if (props.existsPropertyVector<double>(name,
                                               MeshLib::MeshItemType::Node, 1))
        {
            auto const* p = props.getPropertyVector<double>(
                name, MeshLib::MeshItemType::Node, 1);
            auto new_node_vec = new_properties.createNewPropertyVector<double>(
                name, MeshLib::MeshItemType::Node, 1);
            fillNodeProperty(*new_node_vec, *p, node_ids);
            continue;
        }
        if (props.existsPropertyVector<int>(name, MeshLib::MeshItemType::Cell,
                                            1))
        {
            auto const* p = props.getPropertyVector<int>(
                name, MeshLib::MeshItemType::Cell, 1);
            auto new_cell_vec = new_properties.createNewPropertyVector<int>(
                name, MeshLib::MeshItemType::Cell, 1);
            fillElemProperty(*new_cell_vec, *p, elem_ids);
            continue;
        }
        if (props.existsPropertyVector<float>(name, MeshLib::MeshItemType::Cell,
                                              1))
        {
            auto const* p = props.getPropertyVector<float>(
                name, MeshLib::MeshItemType::Cell, 1);
            auto new_cell_vec = new_properties.createNewPropertyVector<float>(
                name, MeshLib::MeshItemType::Cell, 1);
            fillElemProperty(*new_cell_vec, *p, elem_ids);
            continue;
        }
        if (props.existsPropertyVector<double>(name,
                                               MeshLib::MeshItemType::Cell, 1))
        {
            auto const* p = props.getPropertyVector<double>(
                name, MeshLib::MeshItemType::Cell, 1);
            auto new_cell_vec = new_properties.createNewPropertyVector<double>(
                name, MeshLib::MeshItemType::Cell, 1);
            fillElemProperty(*new_cell_vec, *p, elem_ids);
            continue;
        }
        WARN("PropertyVector {:s} not being converted.", name);
    }
    return new_properties;
}
}  // namespace

namespace MeshToolsLib
{
MeshRevision::MeshRevision(MeshLib::Mesh& mesh) : _mesh(mesh) {}

unsigned MeshRevision::getNumberOfCollapsibleNodes(double const eps) const
{
    std::vector<std::size_t> const id_map = collapseNodeIndices(eps);
    std::size_t const nNodes = id_map.size();
    unsigned count(0);
    for (std::size_t i = 0; i < nNodes; ++i)
    {
        if (i != id_map[i])
        {
            count++;
        }
    }
    return count;
}

MeshLib::Mesh* MeshRevision::simplifyMesh(const std::string& new_mesh_name,
                                          double const eps,
                                          unsigned const min_elem_dim) const
{
    if (this->_mesh.getNumberOfElements() == 0)
    {
        return nullptr;
    }

    std::vector<MeshLib::Element*> const& elements(this->_mesh.getElements());
    auto const node_ids = collapseNodeIndices(eps);
    std::vector<MeshLib::Node*> new_nodes =
        this->constructNewNodesArray(node_ids);
    std::vector<MeshLib::Element*> new_elements;
    std::vector<std::size_t> element_ids;

    for (std::size_t k(0); k < elements.size(); ++k)
    {
        MeshLib::Element const* const elem(elements[k]);
        unsigned const n_unique_nodes(getNumberOfUniqueNodes(elem));
        if (n_unique_nodes == elem->getNumberOfBaseNodes() &&
            elem->getDimension() >= min_elem_dim)
        {
            ElementErrorCode const e = elem->validate();
            if (e[ElementErrorFlag::NonCoplanar])
            {
                std::size_t const n_new_elements(
                    subdivideElement(elem, new_nodes, new_elements));
                if (n_new_elements == 0)
                {
                    ERR("Element {:d} has unknown element type.", k);
                    _mesh.resetNodeIDs();
                    BaseLib::cleanupVectorElements(new_nodes, new_elements);
                    return nullptr;
                }
                element_ids.insert(element_ids.end(), n_new_elements, k);
            }
            else
            {
                new_elements.push_back(MeshLib::copyElement(elem, new_nodes));
                element_ids.push_back(k);
            }
        }
        else if (n_unique_nodes < elem->getNumberOfBaseNodes() &&
                 n_unique_nodes > 1)
        {
            std::size_t const n_new_elements(reduceElement(
                elem, n_unique_nodes, new_nodes, new_elements, min_elem_dim));
            element_ids.insert(element_ids.end(), n_new_elements, k);
        }
        else
        {
            ERR("Something is wrong, more unique nodes than actual nodes");
        }
    }

    auto const& props = _mesh.getProperties();
    MeshLib::Properties const new_properties =
        copyProperties(props, node_ids, element_ids);

    _mesh.resetNodeIDs();
    if (!new_elements.empty())
    {
        return new MeshLib::Mesh(new_mesh_name, new_nodes, new_elements,
                                 true /* compute_element_neighbors */,
                                 new_properties);
    }

    BaseLib::cleanupVectorElements(new_nodes, new_elements);
    return nullptr;
}

std::vector<std::size_t> MeshRevision::collapseNodeIndices(
    double const eps) const
{
    const std::vector<MeshLib::Node*>& nodes(_mesh.getNodes());
    const std::size_t nNodes(_mesh.getNumberOfNodes());
    std::vector<std::size_t> id_map(nNodes);
    const double half_eps(eps / 2.0);
    const double sqr_eps(eps * eps);
    std::iota(id_map.begin(), id_map.end(), 0);

    GeoLib::Grid<MeshLib::Node> const grid(nodes.begin(), nodes.end(), 64);

    for (std::size_t k = 0; k < nNodes; ++k)
    {
        MeshLib::Node const* const node(nodes[k]);
        if (node->getID() != k)
        {
            continue;
        }
        std::vector<std::vector<MeshLib::Node*> const*> const node_vectors(
            grid.getPntVecsOfGridCellsIntersectingCube(*node, half_eps));

        const std::size_t nVectors(node_vectors.size());
        for (std::size_t i = 0; i < nVectors; ++i)
        {
            const std::vector<MeshLib::Node*>& cell_vector(*node_vectors[i]);
            const std::size_t nGridCellNodes(cell_vector.size());
            for (std::size_t j = 0; j < nGridCellNodes; ++j)
            {
                MeshLib::Node const* const test_node(cell_vector[j]);
                // are node indices already identical (i.e. nodes will be
                // collapsed)
                if (id_map[node->getID()] == id_map[test_node->getID()])
                {
                    continue;
                }

                // if test_node has already been collapsed to another node x,
                // ignore it (if the current node would need to be collapsed
                // with x it would already have happened when x was tested)
                if (test_node->getID() != id_map[test_node->getID()])
                {
                    continue;
                }

                // calc distance
                if (MathLib::sqrDist(*node, *test_node) < sqr_eps)
                {
                    id_map[test_node->getID()] = node->getID();
                }
            }
        }
    }
    return id_map;
}

std::vector<MeshLib::Node*> MeshRevision::constructNewNodesArray(
    const std::vector<std::size_t>& id_map) const
{
    const std::vector<MeshLib::Node*>& nodes(_mesh.getNodes());
    const std::size_t nNodes(nodes.size());
    std::vector<MeshLib::Node*> new_nodes;
    new_nodes.reserve(nNodes);
    for (std::size_t k = 0; k < nNodes; ++k)
    {
        // all nodes that have not been collapsed with other nodes are copied
        // into new array
        if (nodes[k]->getID() == id_map[k])
        {
            std::size_t const id(new_nodes.size());
            new_nodes.push_back(new MeshLib::Node(
                (*nodes[k])[0], (*nodes[k])[1], (*nodes[k])[2], id));
            nodes[k]->setID(id);  // the node in the old array gets the index of
                                  // the same node in the new array
        }
        // the other nodes are not copied and get the index of the nodes they
        // will have been collapsed with
        else
        {
            nodes[k]->setID(nodes[id_map[k]]->getID());
        }
    }
    return new_nodes;
}

}  // namespace MeshToolsLib
