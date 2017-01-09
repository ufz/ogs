/**
 * \file   AddLayerToMesh.cpp
 * \date   2016-01-18
 * \brief  Implementation of AddLayerToMesh class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AddLayerToMesh.h"

#include <algorithm>
#include <map>
#include <memory>
#include <numeric>
#include <vector>

#include <logog/include/logog.hpp>

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Elements.h"
#include "MeshLib/MeshSurfaceExtraction.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/MeshEditing/FlipElements.h"

namespace MeshLib
{

/** Extrudes point, line, triangle or quad elements to its higher dimensional
 * versions, i.e. line, quad, prism, hexahedron.
 *
 * @param subsfc_nodes the nodes the elements are based on
 * @param sfc_elem the element of the surface that will be extruded
 * @param sfc_to_subsfc_id_map relation between the surface nodes of the surface
 * element and the ids of the nodes of the subsurface mesh
 * @param subsfc_sfc_id_map mapping of the surface nodes of the current mesh
 * to the surface nodes of the extruded mesh
 *
 * @return extruded element (point -> line, line -> quad, tri -> prism, quad ->
 * hexahedron)
*/
MeshLib::Element* extrudeElement(std::vector<MeshLib::Node*> const& subsfc_nodes,
    MeshLib::Element const& sfc_elem,
    MeshLib::PropertyVector<std::size_t> const& sfc_to_subsfc_id_map,
    std::map<std::size_t, std::size_t> const& subsfc_sfc_id_map)
{
    if (sfc_elem.getDimension() > 2)
        return nullptr;

    const unsigned nElemNodes(sfc_elem.getNumberOfBaseNodes());
    MeshLib::Node** new_nodes = new MeshLib::Node*[2*nElemNodes];

    for (unsigned j=0; j<nElemNodes; ++j)
    {
        std::size_t const subsfc_id(
            sfc_to_subsfc_id_map[sfc_elem.getNode(j)->getID()]);
        new_nodes[j] = subsfc_nodes[subsfc_id];
        std::size_t new_idx = (nElemNodes==2) ? (3-j) : (nElemNodes+j);
        new_nodes[new_idx] = subsfc_nodes[subsfc_sfc_id_map.at(subsfc_id)];
    }

    if (sfc_elem.getGeomType() == MeshLib::MeshElemType::LINE)
        return new MeshLib::Quad(new_nodes);
    if (sfc_elem.getGeomType() == MeshLib::MeshElemType::TRIANGLE)
        return new MeshLib::Prism(new_nodes);
    if (sfc_elem.getGeomType() == MeshLib::MeshElemType::QUAD)
        return new MeshLib::Hex(new_nodes);

    return nullptr;
}

MeshLib::Mesh* addTopLayerToMesh(MeshLib::Mesh const& mesh,
    double thickness,
    std::string const& name)
{
    return addLayerToMesh(mesh, thickness, name, true);
}

MeshLib::Mesh* addBottomLayerToMesh(MeshLib::Mesh const& mesh,
    double thickness,
    std::string const& name)
{
    return addLayerToMesh(mesh, thickness, name, false);
}

MeshLib::Mesh* addLayerToMesh(MeshLib::Mesh const& mesh, double thickness,
    std::string const& name,
    bool on_top)
{
    INFO("Extracting top surface of mesh \"%s\" ... ", mesh.getName().c_str());
    int const flag = (on_top) ? -1 : 1;
    const MathLib::Vector3 dir(0, 0, flag);
    double const angle(90);
    std::unique_ptr<MeshLib::Mesh> sfc_mesh (nullptr);

    std::string const prop_name("OriginalSubsurfaceNodeIDs");

    if (mesh.getDimension() == 3)
        sfc_mesh.reset(MeshLib::MeshSurfaceExtraction::getMeshSurface(
            mesh, dir, angle, prop_name));
    else {
        sfc_mesh = (on_top) ? std::unique_ptr<MeshLib::Mesh>(new MeshLib::Mesh(mesh)) :
                              std::unique_ptr<MeshLib::Mesh>(MeshLib::createFlippedMesh(mesh));
        // add property storing node ids
        auto* const pv =
            sfc_mesh->getProperties().createNewPropertyVector<std::size_t>(
                prop_name, MeshLib::MeshItemType::Node, 1);
        if (pv) {
            pv->resize(sfc_mesh->getNumberOfNodes());
            std::iota(pv->begin(), pv->end(), 0);
        } else {
            ERR("Could not create and initialize property.");
            return nullptr;
        }
    }
    INFO("done.");

    // *** add new surface nodes
    std::vector<MeshLib::Node*> subsfc_nodes =
        MeshLib::copyNodeVector(mesh.getNodes());
    std::vector<MeshLib::Element*> subsfc_elements =
        MeshLib::copyElementVector(mesh.getElements(), subsfc_nodes);

    std::size_t const n_subsfc_nodes(subsfc_nodes.size());

    std::vector<MeshLib::Node*> const& sfc_nodes(sfc_mesh->getNodes());
    std::size_t const n_sfc_nodes(sfc_nodes.size());

    // fetch subsurface node ids PropertyVector
    auto const* const node_id_pv =
        sfc_mesh->getProperties().getPropertyVector<std::size_t>(prop_name);
    if (!node_id_pv) {
        ERR(
            "Need subsurface node ids, but the property \"%s\" is not "
            "available.",
            prop_name.c_str());
        return nullptr;
    }

    // *** copy sfc nodes to subsfc mesh node
    std::map<std::size_t, std::size_t> subsfc_sfc_id_map;
    for (std::size_t k(0); k<n_sfc_nodes; ++k) {
        std::size_t const subsfc_id((*node_id_pv)[k]);
        std::size_t const sfc_id(k+n_subsfc_nodes);
        subsfc_sfc_id_map.insert(std::make_pair(subsfc_id, sfc_id));
        MeshLib::Node const& node(*sfc_nodes[k]);
        subsfc_nodes.push_back(new MeshLib::Node(
            node[0], node[1], node[2] - (flag * thickness), sfc_id));
    }

    // *** insert new layer elements into subsfc_mesh
    std::vector<MeshLib::Element*> const& sfc_elements(sfc_mesh->getElements());
    std::size_t const n_sfc_elements(sfc_elements.size());
    for (std::size_t k(0); k<n_sfc_elements; ++k)
        subsfc_elements.push_back(extrudeElement(
            subsfc_nodes, *sfc_elements[k], *node_id_pv, subsfc_sfc_id_map));

    auto new_mesh = new MeshLib::Mesh(name, subsfc_nodes, subsfc_elements);

    auto const* const opt_materials =
        mesh.getProperties().getPropertyVector<int>("MaterialIDs");
    if (!opt_materials)
    {
        ERR("Could not copy the property \"MaterialIDs\" since the original "
            "mesh does not contain such a property.");
        return new_mesh;
    }

    auto* const new_materials =
        new_mesh->getProperties().createNewPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
    if (!new_materials)
    {
        ERR("Can not set material properties for new layer");
        return new_mesh;
    }

    new_materials->reserve(subsfc_elements.size());
    int new_layer_id(
        *(std::max_element(opt_materials->cbegin(), opt_materials->cend())) +
        1);
    std::copy(opt_materials->cbegin(), opt_materials->cend(),
              std::back_inserter(*new_materials));
    auto const n_new_props(subsfc_elements.size() - mesh.getNumberOfElements());
    std::fill_n(std::back_inserter(*new_materials), n_new_props, new_layer_id);

    return new_mesh;
}

} // namespace MeshLib
