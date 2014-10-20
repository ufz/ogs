/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PROCESS_H_
#define PROCESS_LIB_PROCESS_H_

#include "AssemblerLib/MeshComponentMap.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{

class Process
{
public:
    Process(MeshLib::Mesh const& mesh, std::size_t const numberOfNodeComponents)
        : _mesh(mesh), _numberOfNodeComponents(numberOfNodeComponents)
    { }

    virtual ~Process() = default;

    virtual void initialize() = 0;

    std::vector<std::vector<std::size_t>>
    createDofMap(AssemblerLib::MeshComponentMap& mesh_component_map) const
    {
        // dof_map contains for each element vector of global indices to
        // node/element process variables.
        std::vector<std::vector<std::size_t>> dof_map;

        // XXX Shouldn't the mesh_id (and the elements) be from the MeshSubsets?
        // Otherwise we should create the MeshComponentMap also here using the
        // same meshes and elements.
        std::size_t const mesh_id = _mesh.getID();
        std::vector<MeshLib::Element*> const& elements = _mesh.getElements();

        dof_map.reserve(elements.size());

        // For each element find the global indices for node/element components.
        for (MeshLib::Element* e : elements)
        {
            std::size_t const nnodes = e->getNNodes();
            std::vector<MeshLib::Location> vec_items;
            vec_items.reserve(nnodes);
            for (std::size_t j = 0; j < nnodes; j++)
                vec_items.emplace_back(
                    mesh_id,
                    MeshLib::MeshItemType::Node,
                    e->getNode(j)->getID());

            dof_map.push_back(mesh_component_map.getGlobalIndices<
                AssemblerLib::ComponentOrder::BY_COMPONENT>(vec_items));
        }

        return dof_map;
    }

protected:
    MeshLib::Mesh const& _mesh;
    std::size_t const _numberOfNodeComponents;
};

}   // namespace ProcessLib

//
// Include all known processes here.
//
#include "GroundwaterFlowProcess.h"

#endif  // PROCESS_LIB_PROCESS_H_
