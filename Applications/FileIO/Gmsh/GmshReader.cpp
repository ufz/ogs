/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "GmshReader.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <map>
#include <type_traits>
#include <vector>

#include "BaseLib/FileTools.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshEditing/ElementValueModification.h"

namespace FileIO
{
namespace GMSH
{
bool isGMSHMeshFile(const std::string& fname)
{
    std::ifstream input(fname.c_str());

    if (!input)
    {
        ERR("isGMSHMeshFile(): Could not open file {:s}.", fname);
        return false;
    }

    std::string header_first_line;
    input >> header_first_line;
    if (header_first_line.find("$MeshFormat") != std::string::npos)
    {
        // read version
        std::string version;
        std::getline(input, version);
        std::getline(input, version);
        INFO("isGMSHMeshFile(): Found GMSH mesh file version: {:s}.", version);
        input.close();
        return true;
    }

    return false;
}

void readNodeIDs(std::ifstream& in, unsigned n_nodes,
                 std::vector<unsigned>& node_ids,
                 std::map<unsigned, unsigned> const& id_map)
{
    unsigned idx;
    for (unsigned i = 0; i < n_nodes; i++)
    {
        in >> idx;
        node_ids.push_back(id_map.at(idx));
    }
}

template <typename ElementType>
std::pair<MeshLib::Element*, int> createElement(
    std::ifstream& in, std::vector<MeshLib::Node*> const& nodes,
    int const mat_id, std::map<unsigned, unsigned> const& id_map)
{
    std::vector<unsigned> node_ids;
    readNodeIDs(in, ElementType::n_all_nodes, node_ids, id_map);

    std::array<MeshLib::Node*, ElementType::n_all_nodes> element_nodes;

    std::transform(begin(node_ids), end(node_ids), begin(element_nodes),
                   [&nodes](auto const id) { return nodes[id]; });

    return std::make_pair(new ElementType(element_nodes), mat_id);
}

template <>
std::pair<MeshLib::Element*, int> createElement<MeshLib::Tri>(
    std::ifstream& in, std::vector<MeshLib::Node*> const& nodes,
    int const mat_id, std::map<unsigned, unsigned> const& id_map)
{
    std::vector<unsigned> node_ids;
    readNodeIDs(in, 3, node_ids, id_map);

    std::array<MeshLib::Node*, 3> element_nodes;

    std::transform(std::rbegin(node_ids), std::rend(node_ids),
                   begin(element_nodes),
                   [&nodes](auto const id) { return nodes[id]; });

    return std::make_pair(new MeshLib::Tri(element_nodes), mat_id);
}

template <>
std::pair<MeshLib::Element*, int> createElement<MeshLib::Tet10>(
    std::ifstream& in, std::vector<MeshLib::Node*> const& nodes,
    int const mat_id, std::map<unsigned, unsigned> const& id_map)
{
    std::vector<unsigned> node_ids;
    readNodeIDs(in, MeshLib::Tet10::n_all_nodes, node_ids, id_map);

    std::swap(node_ids[8], node_ids[9]);

    std::array<MeshLib::Node*, MeshLib::Tet10::n_all_nodes> element_nodes;

    std::transform(begin(node_ids), end(node_ids), begin(element_nodes),
                   [&nodes](auto const id) { return nodes[id]; });

    return std::make_pair(new MeshLib::Tet10(element_nodes), mat_id);
}

template <>
std::pair<MeshLib::Element*, int> createElement<MeshLib::Hex20>(
    std::ifstream& in, std::vector<MeshLib::Node*> const& nodes,
    int const mat_id, std::map<unsigned, unsigned> const& id_map)
{
    std::vector<unsigned> node_ids;
    readNodeIDs(in, MeshLib::Hex20::n_all_nodes, node_ids, id_map);

    std::array<MeshLib::Node*, MeshLib::Hex20::n_all_nodes> element_nodes;

    constexpr std::array node_order = {0,  1, 2,  3,  4,  5,  6,  7,  8,  11,
                                       13, 9, 16, 18, 19, 17, 10, 12, 14, 15};

    std::transform(begin(node_order), end(node_order), begin(element_nodes),
                   [&node_ids, &nodes](auto const id)
                   { return nodes[node_ids[id]]; });

    return std::make_pair(new MeshLib::Hex20(element_nodes), mat_id);
}

template <>
std::pair<MeshLib::Element*, int> createElement<MeshLib::Prism15>(
    std::ifstream& in, std::vector<MeshLib::Node*> const& nodes,
    int const mat_id, std::map<unsigned, unsigned> const& id_map)
{
    std::vector<unsigned> node_ids;
    readNodeIDs(in, MeshLib::Prism15::n_all_nodes, node_ids, id_map);

    std::array<MeshLib::Node*, MeshLib::Prism15::n_all_nodes> element_nodes;

    constexpr std::array node_order = {0, 1,  2,  3,  4, 5,  6, 9,
                                       7, 12, 14, 13, 8, 10, 11};

    std::transform(begin(node_order), end(node_order), begin(element_nodes),
                   [&node_ids, &nodes](auto const id)
                   { return nodes[node_ids[id]]; });

    return std::make_pair(new MeshLib::Prism15(element_nodes), mat_id);
}

template <>
std::pair<MeshLib::Element*, int> createElement<MeshLib::Pyramid13>(
    std::ifstream& in, std::vector<MeshLib::Node*> const& nodes,
    int const mat_id, std::map<unsigned, unsigned> const& id_map)
{
    std::vector<unsigned> node_ids;
    readNodeIDs(in, MeshLib::Pyramid13::n_all_nodes, node_ids, id_map);
    std::array<MeshLib::Node*, MeshLib::Pyramid13::n_all_nodes> element_nodes;

    constexpr std::array node_order = {0,  1, 2, 3, 4,  5, 8,
                                       10, 6, 7, 9, 11, 12};

    std::transform(begin(node_order), end(node_order), begin(element_nodes),
                   [&node_ids, &nodes](auto const id)
                   { return nodes[node_ids[id]]; });

    return std::make_pair(new MeshLib::Pyramid13(element_nodes), mat_id);
}

std::pair<MeshLib::Element*, int> readElement(
    std::ifstream& in, std::vector<MeshLib::Node*> const& nodes,
    std::map<unsigned, unsigned> const& id_map)
{
    unsigned idx;
    unsigned type;
    unsigned n_tags;
    unsigned dummy;
    int mat_id;

    // element format is structured like this:
    // element-id element-type n-tags physical-entity elementary entity node-ids
    in >> idx >> type >> n_tags >> dummy >> mat_id;

    switch (type)
    {
        case 1:
        {
            return createElement<MeshLib::Line>(in, nodes, mat_id, id_map);
        }
        case 2:
        {
            return createElement<MeshLib::Tri>(in, nodes, mat_id, id_map);
        }
        case 3:
        {
            return createElement<MeshLib::Quad>(in, nodes, mat_id, id_map);
        }
        case 4:
        {
            return createElement<MeshLib::Tet>(in, nodes, mat_id, id_map);
        }
        case 5:
        {
            return createElement<MeshLib::Hex>(in, nodes, mat_id, id_map);
        }
        case 6:
        {
            return createElement<MeshLib::Prism>(in, nodes, mat_id, id_map);
        }
        case 7:
        {
            return createElement<MeshLib::Pyramid>(in, nodes, mat_id, id_map);
        }
        case 8:  // 3-node second order line.
        {
            return createElement<MeshLib::Line3>(in, nodes, mat_id, id_map);
        }
        case 9:  // 6-node second order triangle.
        {
            return createElement<MeshLib::Tri6>(in, nodes, mat_id, id_map);
        }
        case 10:  // 9-node second order quadrangle.
        {
            return createElement<MeshLib::Quad9>(in, nodes, mat_id, id_map);
        }
        case 11:  // 10-node second order tetrahedron.
        {
            return createElement<MeshLib::Tet10>(in, nodes, mat_id, id_map);
        }
        case 16:  // 8-node second order quadrangle.
        {
            return createElement<MeshLib::Quad8>(in, nodes, mat_id, id_map);
        }
        case 17:  // 20-node second order hexahedron.
        {
            return createElement<MeshLib::Hex20>(in, nodes, mat_id, id_map);
        }
        case 18:  // 15-node second order prism.
        {
            return createElement<MeshLib::Prism15>(in, nodes, mat_id, id_map);
        }
        case 19:  // 13-node second order pyramid.
        {
            return createElement<MeshLib::Pyramid13>(in, nodes, mat_id, id_map);
        }
        case 15:
            in >> dummy;  // skip rest of line
            break;
        default:
            WARN("readGMSHMesh(): Unknown element type {:d}.", type);
            break;
    }
    return std::make_pair(nullptr, -1);
}

MeshLib::Mesh* readGMSHMesh(std::string const& fname)
{
    std::string line;
    std::ifstream in(fname.c_str(), std::ios::in);
    if (!in.is_open())
    {
        WARN("readGMSHMesh() - Could not open file {:s}.", fname);
        return nullptr;
    }

    std::getline(in, line);  // $MeshFormat keyword
    if (line.find("$MeshFormat") == std::string::npos)
    {
        in.close();
        WARN("No GMSH file format recognized.");
        return nullptr;
    }

    std::getline(in, line);  // version-number file-type data-size
    if (line.substr(0, 3) != "2.2")
    {
        WARN("Wrong gmsh file format version '{:s}'.", line.substr(0, 3));
        return nullptr;
    }

    if (line[4] != '0')
    {
        WARN("Currently reading gmsh binary file type is not supported.");
        return nullptr;
    }
    std::getline(in, line);  //$EndMeshFormat

    std::vector<MeshLib::Node*> nodes;
    std::vector<MeshLib::Element*> elements;
    std::vector<int> materials;
    std::map<unsigned, unsigned> id_map;
    while (line.find("$EndElements") == std::string::npos)
    {
        // Node data
        std::getline(in, line);  //$Nodes Keywords
        if (line.find("$Nodes") != std::string::npos)
        {
            std::size_t n_nodes(0);
            long id;
            double x;
            double y;
            double z;
            in >> n_nodes >> std::ws;
            nodes.resize(n_nodes);
            for (std::size_t i = 0; i < n_nodes; i++)
            {
                in >> id >> x >> y >> z >> std::ws;
                id_map.insert(std::map<unsigned, unsigned>::value_type(id, i));
                nodes[i] = new MeshLib::Node(x, y, z, id);
            }
            std::getline(in, line);  // End Node keyword $EndNodes
        }

        // Element data
        if (line.find("$Elements") != std::string::npos)
        {
            std::size_t n_elements(0);
            if (!(in >> n_elements >> std::ws))
            {  // number-of-elements
                ERR("Read GMSH mesh does not contain any elements");
            }
            elements.reserve(n_elements);
            materials.reserve(n_elements);
            for (std::size_t i = 0; i < n_elements; i++)
            {
                MeshLib::Element* elem(nullptr);
                int mat_id(0);
                std::tie(elem, mat_id) = readElement(in, nodes, id_map);

                if (elem)
                {
                    elements.push_back(elem);
                    materials.push_back(mat_id);
                }
            }
            std::getline(in, line);  // END keyword
        }

        if (line.find("PhysicalNames") != std::string::npos)
        {
            std::size_t n_lines(0);
            in >> n_lines >> std::ws;  // number-of-lines
            for (std::size_t i = 0; i < n_lines; i++)
            {
                std::getline(in, line);
            }
            std::getline(in, line);  // END keyword
        }
    }
    in.close();
    if (elements.empty())
    {
        for (auto& node : nodes)
        {
            delete node;
        }
        return nullptr;
    }

    MeshLib::Mesh* mesh(new MeshLib::Mesh(
        BaseLib::extractBaseNameWithoutExtension(fname), nodes, elements,
        true /* compute_element_neighbors */));

    auto* const material_ids =
        mesh->getProperties().createNewPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
    if (!material_ids)
    {
        WARN("Could not create PropertyVector for MaterialIDs in Mesh.");
    }
    else
    {
        material_ids->insert(material_ids->end(), materials.cbegin(),
                             materials.cend());
    }

    MeshToolsLib::ElementValueModification::condense(*mesh);

    INFO("\t... finished.");
    INFO("Nr. Nodes: {:d}.", nodes.size());
    INFO("Nr. Elements: {:d}.", elements.size());

    return mesh;
}

}  // end namespace GMSH
}  // end namespace FileIO
