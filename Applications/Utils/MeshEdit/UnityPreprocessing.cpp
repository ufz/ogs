/*
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <string>
#include <vector>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"



bool containsCellVecs(MeshLib::Mesh const& mesh)
{
    MeshLib::Properties const& props (mesh.getProperties());
    std::vector<std::string> const& vec_names(props.getPropertyVectorNames(MeshLib::MeshItemType::Cell));
    if (vec_names.empty())
        return false;
    return true;
}

template <class T>
MeshLib::Element* createElement(MeshLib::Element const& e,
    std::vector<MeshLib::Node*> &nodes,
    std::vector<std::vector<std::size_t>> &node_map)
{
    std::size_t const n_nodes = e.getNumberOfNodes();
    MeshLib::Node** new_nodes = new MeshLib::Node*[n_nodes];
    for (unsigned i = 0; i<n_nodes; ++i)
    {
        new_nodes[i] = new MeshLib::Node(e.getNode(i)->getCoords(), nodes.size());
        node_map[e.getNodeIndex(i)].push_back(nodes.size());
        nodes.push_back(new_nodes[i]);
    }
    return new T(new_nodes);
}

template <class T>
bool fillPropVec(MeshLib::Properties const& props,
    std::string const& name,
    MeshLib::Properties &new_props,
    std::vector<MeshLib::Element*> const& elems,
    std::vector<std::vector<std::size_t>> const& node_map,
    std::size_t total_nodes)
{
    if (!props.existsPropertyVector<T>(name))
        return false;

    MeshLib::PropertyVector<T> const*const vec = props.getPropertyVector<T>(name);
    if (vec->getNumberOfComponents() != 1)
    {
        INFO("Ignoring array \"%s\" (more than one component).", name.c_str());
        return false;
    }

    MeshLib::PropertyVector<T>* new_vec =
        new_props.createNewPropertyVector<T>(vec->getPropertyName(), MeshLib::MeshItemType::Node, 1);
    new_vec->resize(total_nodes);
    if (vec->getMeshItemType() == MeshLib::MeshItemType::Node)
    {
        INFO("Migrating node array \"%s\" to new mesh structure...", name.c_str());
        std::size_t const n_nodes (node_map.size());
        for (std::size_t i = 0; i<n_nodes; ++i)
        {
            std::size_t const n_nodes = node_map[i].size();
            for (std::size_t j = 0; j<n_nodes; ++j)
                (*new_vec)[node_map[i][j]] = (*vec)[i];
        }
    }
    else if (vec->getMeshItemType() == MeshLib::MeshItemType::Cell)
    {
        INFO("Transforming cell array \"%s\" into node array...", name.c_str());
        std::size_t const n_elems (vec->size());
        for (std::size_t i = 0; i<n_elems; ++i)
        {
            std::size_t const n_nodes = elems[i]->getNumberOfNodes();
            for (std::size_t j = 0; j<n_nodes; ++j)
                (*new_vec)[elems[i]->getNodeIndex(j)] = (*vec)[i];
        }
    }
    return true;
}

MeshLib::Properties constructProperties(MeshLib::Properties const& props,
    std::vector<MeshLib::Element*> const& elems,
    std::vector<std::vector<std::size_t>> const& node_map,
    std::size_t n_nodes)
{
    std::vector<std::string> const& names = props.getPropertyVectorNames();
    MeshLib::Properties new_props;
    for (std::string const& name : names)
    {
        if (fillPropVec<int>(props, name, new_props, elems, node_map, n_nodes))
            continue;
        if (fillPropVec<double>(props, name, new_props, elems, node_map, n_nodes))
            continue;
        if (fillPropVec<long>(props, name, new_props, elems, node_map, n_nodes))
            continue;
    }
    return new_props;
}

MeshLib::Mesh* constructMesh(MeshLib::Mesh const& mesh)
{
    INFO("Splitting nodes...");
    std::vector<MeshLib::Element*> const& elems = mesh.getElements();
    std::vector<MeshLib::Node*> new_nodes;
    std::vector<MeshLib::Element*> new_elems;
    std::vector<std::vector<std::size_t>> node_map;
    node_map.resize(mesh.getNumberOfNodes());
    for (MeshLib::Element* elem : elems)
    {
        if (elem->getGeomType() == MeshLib::MeshElemType::LINE)
            new_elems.push_back(createElement<MeshLib::Line>(*elem, new_nodes, node_map));
        else if (elem->getGeomType() == MeshLib::MeshElemType::TRIANGLE)
            new_elems.push_back(createElement<MeshLib::Tri>(*elem, new_nodes, node_map));
        else if (elem->getGeomType() == MeshLib::MeshElemType::QUAD)
            new_elems.push_back(createElement<MeshLib::Quad>(*elem, new_nodes, node_map));
        else if (elem->getGeomType() == MeshLib::MeshElemType::TETRAHEDRON)
            new_elems.push_back(createElement<MeshLib::Tet>(*elem, new_nodes, node_map));
        else if (elem->getGeomType() == MeshLib::MeshElemType::HEXAHEDRON)
            new_elems.push_back(createElement<MeshLib::Hex>(*elem, new_nodes, node_map));
        else if (elem->getGeomType() == MeshLib::MeshElemType::PYRAMID)
            new_elems.push_back(createElement<MeshLib::Pyramid>(*elem, new_nodes, node_map));
        else if (elem->getGeomType() == MeshLib::MeshElemType::PRISM)
            new_elems.push_back(createElement<MeshLib::Prism>(*elem, new_nodes, node_map));
        else
        {
            ERR("Error: Unknown element type.");
            return nullptr;
        }
    }

    MeshLib::Properties new_props =
        constructProperties(mesh.getProperties(), new_elems, node_map, new_nodes.size());
    return new MeshLib::Mesh("Unity conform mesh", new_nodes, new_elems, new_props);
}

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Prepares OGS-meshes for use in Unity"
        "",
        ' ',
        "0.1");

    TCLAP::ValueArg<std::string> mesh_arg("i", "input",
        "the file containing the original OGS mesh", true,
        "", "input file name");
    cmd.add(mesh_arg);

    TCLAP::ValueArg<std::string> mesh_out_arg("o", "output",
        "the file name the result will be written to", true,
        "", "output file name");
    cmd.add(mesh_out_arg);

    cmd.parse(argc, argv);

    INFO("Reading mesh \"%s\" ... ", mesh_arg.getValue().c_str());
    std::unique_ptr<MeshLib::Mesh> mesh {MeshLib::IO::readMeshFromFile(mesh_arg.getValue())};
    if (!mesh)
        return EXIT_FAILURE;
    INFO("done.\n");

    INFO("Checking for line elements...");
    std::array<unsigned, 7> const& n_element_types =
        MeshLib::MeshInformation::getNumberOfElementTypes(*mesh);
    std::unique_ptr<MeshLib::Mesh> result;
    if (n_element_types[0] == 0)
    {
        INFO ("No line elements found.\n");
        result = std::move(mesh);
    }
    else if (n_element_types[0] == mesh->getNumberOfElements())
    {
        INFO ("Keeping line mesh.\n");
        result = std::move(mesh);
    }
    else
    {
        MeshLib::ElementSearch searcher(*mesh);
        std::size_t const n_rem_elems = searcher.searchByElementType(MeshLib::MeshElemType::LINE);
        result.reset(MeshLib::removeElements(*mesh, searcher.getSearchedElementIDs(), "temp mesh"));
        INFO ("%d line elements found and removed.\n", n_rem_elems);
    }

    INFO("Checking for cell-arrays...");
    if (containsCellVecs(*result))
        result.reset(constructMesh(*result));
    else
        INFO("No cell arrays found, keeping mesh structure.\n");

    INFO("Writing mesh \"%s\" ... ", mesh_out_arg.getValue().c_str());
    MeshLib::IO::VtuInterface writer(result.get(), vtkXMLWriter::Ascii, false);
    writer.writeToFile(mesh_out_arg.getValue());
    INFO("done.");

    return EXIT_SUCCESS;
}
