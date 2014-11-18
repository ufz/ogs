/**
 * \file LocalToGlobalIndexMap.cpp
 * \author Norihiro Watanabe
 * \author Wenqing Wang
 * \date   2013-04-16, 2014-11-14
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LocalToGlobalIndexMap.h"

#include "logog/include/logog.hpp"

#include "AssemblerLib/MeshComponentMap.h"
#include "MeshLib/MeshSubsets.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"
#endif

namespace AssemblerLib
{

#ifndef USE_PETSC
LocalToGlobalIndexMap::LocalToGlobalIndexMap(
    std::vector<MeshLib::MeshSubsets*> const& mesh_subsets,
    AssemblerLib::ComponentOrder const order, const bool is_linear_element)
    : _mesh_subsets(mesh_subsets), _mesh_component_map(_mesh_subsets, order)
{
    // For all MeshSubsets and each of their MeshSubset's and each element
    // of that MeshSubset save a line of global indices.
    for (MeshLib::MeshSubsets const* const mss : _mesh_subsets)
    {
        for (MeshLib::MeshSubset const* const ms : *mss)
        {
            std::size_t const mesh_id = ms->getMeshID();

            // For each element find the global indices for node/element
            // components.
            for (auto e = ms->elementsBegin();
                    e != ms->elementsEnd(); ++e)
            {
                std::vector<MeshLib::Location> vec_items;
                std::size_t nnodes = (*e)->getNNodes();
                if(is_is_linear_element)
                    nnodes = (*e)->getNBaseNodes();
                vec_items.reserve(nnodes);

                for (std::size_t n = 0; n < nnodes; n++)
                {
                    vec_items.emplace_back(
                        mesh_id,
                        MeshLib::MeshItemType::Node,
                        (*e)->getNode(n)->getID());
                }

                // Save a line of indices for the current element.
                switch (order)
                {
                    case AssemblerLib::ComponentOrder::BY_LOCATION:
                        _rows.push_back(_mesh_component_map.getGlobalIndices<AssemblerLib::ComponentOrder::BY_LOCATION>(vec_items));
                        break;
                    case AssemblerLib::ComponentOrder::BY_COMPONENT:
                        _rows.push_back(_mesh_component_map.getGlobalIndices<AssemblerLib::ComponentOrder::BY_COMPONENT>(vec_items));
                        break;
                }
            }
        }
    }    
}

#else   //#if defined(USE_PETSC)
LocalToGlobalIndexMap::LocalToGlobalIndexMap(
    std::vector<MeshLib::MeshSubsets*> const& mesh_subsets,
    AssemblerLib::ComponentOrder const order, const bool is_linear_element)
    : _mesh_subsets(mesh_subsets), _mesh_component_map(_mesh_subsets, order)
{
    // For all MeshSubsets and each of their MeshSubset's and each element
    // of that MeshSubset save a line of global indices.
    for (MeshLib::MeshSubsets const* const mss : _mesh_subsets)
    {
        for (MeshLib::MeshSubset const* const ms : *mss)
        {
            std::size_t const mesh_id = ms->getMeshID();

            const MeshLib::NodePartitionedMesh &mesh 
                    = dynamic_cast<const MeshLib::NodePartitionedMesh&>(ms->getMesh());
                        
            // For each element find the global indices for node/element
            // components.
            for (auto e = ms->elementsBegin();
                    e != ms->elementsEnd(); ++e)
            {
                std::vector<MeshLib::Location> vec_items;
                std::size_t nnodes = (*e)->getNNodes();
                if(is_linear_element)
                    nnodes = (*e)->getNBaseNodes();
                vec_items.reserve(nnodes);
                std::vector<bool> ghost_node_flag(nnodes+1);

                ghost_node_flag[0] = false; // Indicator for ghost element		
                for (std::size_t n = 0; n < nnodes; n++)
                {
                    const size_t node_id = (*e)->getNode(n)->getID(); 
                    const size_t global_node_id = mesh.getGlobalNodeID(node_id);                     					
                    vec_items.emplace_back(
                        mesh_id,
                        MeshLib::MeshItemType::Node, global_node_id);
                        
                    if( mesh.isGhostNode(node_id) )
                    {                   
                       ghost_node_flag[n+1] = true;
                       ghost_node_flag[0] = true;                       
                    }                        
                    else
                       ghost_node_flag[n+1] = false;                            
                }
                
                _element_ghost_node_flags.push_back(ghost_node_flag);

                // Save a line of indices for the current element.
                switch (order)
                {
                    case AssemblerLib::ComponentOrder::BY_LOCATION:
                        {                    
                            _rows.push_back(_mesh_component_map.getGlobalIndices<AssemblerLib::ComponentOrder::BY_LOCATION>(vec_items));
                        }
                        break;
                    case AssemblerLib::ComponentOrder::BY_COMPONENT:
                        {                     
                            _rows.push_back(_mesh_component_map.getGlobalIndices<AssemblerLib::ComponentOrder::BY_COMPONENT>(vec_items));
                        } 
                        break;
                }
            }
        }
    }    
} 
#endif  //#ifndef USE_PETSC

std::size_t
LocalToGlobalIndexMap::dofSize() const
{
    return _mesh_component_map.size();
}

std::size_t
LocalToGlobalIndexMap::size() const
{
    return _rows.size();
}

LocalToGlobalIndexMap::RowColumnIndices
LocalToGlobalIndexMap::operator[](std::size_t const mesh_item_id) const
{
    return RowColumnIndices(_rows[mesh_item_id], _columns[mesh_item_id]);
}

LocalToGlobalIndexMap::LineIndex
LocalToGlobalIndexMap::rowIndices(std::size_t const mesh_item_id) const
{
    return _rows[mesh_item_id];
}

LocalToGlobalIndexMap::LineIndex
LocalToGlobalIndexMap::columnIndices(std::size_t const mesh_item_id) const
{
    return _columns[mesh_item_id];
}

}   // namespace AssemblerLib
