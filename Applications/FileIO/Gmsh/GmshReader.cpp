/**
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "GmshReader.h"

#include <fstream>
#include <map>
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
#include "MeshLib/MeshEditing/ElementValueModification.h"
#include "MeshLib/Node.h"

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
        getline(input, version);
        getline(input, version);
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

std::pair<MeshLib::Element*, int> readElement(
    std::ifstream& in, std::vector<MeshLib::Node*> const& nodes,
    std::map<unsigned, unsigned> const& id_map)
{
    unsigned idx;
    unsigned type;
    unsigned n_tags;
    unsigned dummy;
    int mat_id;
    std::vector<unsigned> node_ids;

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

    getline(in, line);  // $MeshFormat keyword
    if (line.find("$MeshFormat") == std::string::npos)
    {
        in.close();
        WARN("No GMSH file format recognized.");
        return nullptr;
    }

    getline(in, line);  // version-number file-type data-size
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
    getline(in, line);  //$EndMeshFormat

    std::vector<MeshLib::Node*> nodes;
    std::vector<MeshLib::Element*> elements;
    std::vector<int> materials;
    std::map<unsigned, unsigned> id_map;
    while (line.find("$EndElements") == std::string::npos)
    {
        // Node data
        getline(in, line);  //$Nodes Keywords
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
            getline(in, line);  // End Node keyword $EndNodes
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
            getline(in, line);  // END keyword
        }

        if (line.find("PhysicalNames") != std::string::npos)
        {
            std::size_t n_lines(0);
            in >> n_lines >> std::ws;  // number-of-lines
            for (std::size_t i = 0; i < n_lines; i++)
            {
                getline(in, line);
            }
            getline(in, line);  // END keyword
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
        BaseLib::extractBaseNameWithoutExtension(fname), nodes, elements));

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

    MeshLib::ElementValueModification::condense(*mesh);

    INFO("\t... finished.");
    INFO("Nr. Nodes: {:d}.", nodes.size());
    INFO("Nr. Elements: {:d}.", elements.size());

    return mesh;
}

}  // end namespace GMSH
}  // end namespace FileIO
