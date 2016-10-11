/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "GmshReader.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Pyramid.h"

#include "MeshLib/Node.h"
#include "MeshLib/Mesh.h"

#include "MeshLib/MeshEditing/ElementValueModification.h"

#include "BaseLib/FileTools.h"

#include <fstream>
#include <map>
#include <vector>

namespace FileIO
{
namespace GMSH
{

bool isGMSHMeshFile(const std::string& fname)
{
    std::ifstream input(fname.c_str());

    if (!input) {
        ERR("isGMSHMeshFile(): Could not open file %s.", fname.c_str());
        return false;
    }

    std::string header_first_line;
    input >> header_first_line;
    if (header_first_line.find("$MeshFormat") != std::string::npos) {
        // read version
        std::string version;
        getline(input, version);
        getline(input, version);
        INFO("isGMSHMeshFile(): Found GMSH mesh file version: %s.",
             version.c_str());
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

std::pair<MeshLib::Element*, int> readElement(
    std::ifstream& in, std::vector<MeshLib::Node*> const& nodes,
    std::map<unsigned, unsigned> const& id_map)
{
    unsigned idx, type, n_tags, dummy;
    int mat_id;
    std::vector<unsigned> node_ids;
    std::vector<MeshLib::Node*> elem_nodes;
    in >> idx >> type >> n_tags >> mat_id >> dummy;

    // skip tags
    for (std::size_t j = 2; j < n_tags; j++)
        in >> dummy;

    switch (type)
    {
    case 1: {
        readNodeIDs(in, 2, node_ids, id_map);
        // edge_nodes array will be deleted from Line object
        auto edge_nodes = new MeshLib::Node*[2];
        edge_nodes[0] = nodes[node_ids[0]];
        edge_nodes[1] = nodes[node_ids[1]];
        return std::make_pair(new MeshLib::Line(edge_nodes), mat_id);
    }
    case 2: {
        readNodeIDs(in, 3, node_ids, id_map);
        auto tri_nodes = new MeshLib::Node*[3];
        tri_nodes[0] = nodes[node_ids[2]];
        tri_nodes[1] = nodes[node_ids[1]];
        tri_nodes[2] = nodes[node_ids[0]];
        return std::make_pair(new MeshLib::Tri(tri_nodes), mat_id);
    }
    case 3: {
        readNodeIDs(in, 4, node_ids, id_map);
        auto quad_nodes = new MeshLib::Node*[4];
        for (unsigned k(0); k < 4; k++)
            quad_nodes[k] = nodes[node_ids[k]];
        return std::make_pair(new MeshLib::Quad(quad_nodes), mat_id);
    }
    case 4: {
        readNodeIDs(in, 4, node_ids, id_map);
        auto tet_nodes = new MeshLib::Node*[5];
        for (unsigned k(0); k < 4; k++)
            tet_nodes[k] = nodes[node_ids[k]];
        return std::make_pair(new MeshLib::Tet(tet_nodes), mat_id);
    }
    case 5: {
        readNodeIDs(in, 8, node_ids, id_map);
        auto hex_nodes = new MeshLib::Node*[8];
        for (unsigned k(0); k < 8; k++)
            hex_nodes[k] = nodes[node_ids[k]];
        return std::make_pair(new MeshLib::Hex(hex_nodes), mat_id);
    }
    case 6: {
        readNodeIDs(in, 6, node_ids, id_map);
        auto prism_nodes = new MeshLib::Node*[6];
        for (unsigned k(0); k < 6; k++)
            prism_nodes[k] = nodes[node_ids[k]];
        return std::make_pair(new MeshLib::Prism(prism_nodes), mat_id);
    }
    case 7: {
        readNodeIDs(in, 5, node_ids, id_map);
        auto pyramid_nodes = new MeshLib::Node*[5];
        for (unsigned k(0); k < 5; k++)
            pyramid_nodes[k] = nodes[node_ids[k]];
        return std::make_pair(new MeshLib::Pyramid(pyramid_nodes), mat_id);
    }
    case 15:
        in >> dummy; // skip rest of line
        break;
    default:
        WARN("readGMSHMesh(): Unknown element type %d.", type);
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
        WARN ("readGMSHMesh() - Could not open file %s.", fname.c_str());
        return nullptr;
    }

    getline(in, line); // $MeshFormat keyword
    if (line.find("$MeshFormat") == std::string::npos)
    {
        in.close();
        WARN ("No GMSH file format recognized.");
        return nullptr;
    }

    getline(in, line); // version-number file-type data-size
    if (line.substr(0,3).compare("2.2") != 0) {
        WARN("Wrong gmsh file format version.");
        return nullptr;
    }

    if (line.substr(4,1).compare("0") != 0) {
        WARN("Currently reading gmsh binary file type is not supported.");
        return nullptr;
    }
    getline(in, line); //$EndMeshFormat

    std::vector<MeshLib::Node*> nodes;
    std::vector<MeshLib::Element*> elements;
    std::vector<int> materials;
    std::map<unsigned, unsigned> id_map;
    while (line.find("$EndElements") == std::string::npos)
    {
        // Node data
        getline(in, line); //$Nodes Keywords
        if (line.find("$Nodes") != std::string::npos)
        {
            std::size_t n_nodes(0);
            long id;
            double x, y, z;
            in >> n_nodes >> std::ws;
            nodes.resize(n_nodes);
            for (std::size_t i = 0; i < n_nodes; i++) {
                in >> id >> x >> y >> z >> std::ws;
                id_map.insert(std::map<unsigned, unsigned>::value_type(id, i));
                nodes[i] = new MeshLib::Node(x,y,z,id);
            }
            getline(in, line); // End Node keyword $EndNodes
        }

        // Element data
        if (line.find("$Elements") != std::string::npos)
        {
            std::size_t n_elements(0);
            if (! (in >> n_elements >> std::ws)) { // number-of-elements
                ERR("Read GMSH mesh does not contain any elements");
            }
            elements.reserve(n_elements);
            materials.reserve(n_elements);
            for (std::size_t i = 0; i < n_elements; i++)
            {
                MeshLib::Element* elem(nullptr);
                int mat_id(0);
                std::tie(elem, mat_id) = readElement(in, nodes, id_map);

                if (elem) {
                    elements.push_back(elem);
                    materials.push_back(mat_id);
                }
            }
            getline(in, line); // END keyword
        }

        if (line.find("PhysicalNames") != std::string::npos)
        {
            std::size_t n_lines(0);
            in >> n_lines >> std::ws; // number-of-lines
            for (std::size_t i = 0; i < n_lines; i++)
                getline(in, line);
            getline(in, line); // END keyword
        }
    }
    in.close();
    if (elements.empty()) {
        for (auto it(nodes.begin()); it != nodes.end(); ++it) {
            delete *it;
        }
        return nullptr;
    }

    MeshLib::Mesh * mesh(new MeshLib::Mesh(
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
    INFO("Nr. Nodes: %d.", nodes.size());
    INFO("Nr. Elements: %d.", elements.size());

    return mesh;
}

} // end namespace GMSH
} // end namespace FileIO
