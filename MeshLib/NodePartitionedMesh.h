/*!
  \file NodePartitionedMesh.h
  \author Wenqing Wang
  \date   2014.06
  \brief  Definition of mesh class for partitioned mesh (by node) for parallel computing within the
          framework of domain decomposition (DDC).

  \copyright
  Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#ifndef NODE_PARTITIONED_MESH_H_
#define NODE_PARTITIONED_MESH_H_

#include <vector>
#include <string>

#include "Mesh.h"

namespace FileIO
{
class readNodePartitionedMesh;
};

namespace MeshLib
{
class Node;
class Element;

/// A subdomain mesh.
class NodePartitionedMesh : public Mesh
{
    public:
        /*!
            \brief Constructor
            \param name          Name assigned to the mesh.
            \param nodes         Vector for nodes, which storage looks like: 
                                 ||--active base nodes--|--ghost base nodes--|
                                  --active extra nodes--|--ghost extra nodes--||
                                 (extra nodes: nodes for high order interpolations)  
            \param glb_node_ids  Global IDs of nodes of a partition.
            \param elements      Vector for elements. Ghost elements are stored
                                 after regular (non-ghost) elements.
            \param g_elem_data   Element wise local node IDs of active nodes.
            \param n_nghost_elem Number of non-ghost elements, or the start ID of
                                 the entry of ghost element in the element vector.
            \param nnodes_global Number of nodes of the global mesh.
                                 0: with linear elements
                                 1: with quadratic elemens.
            \param nnodes_active Number of active nodes of the partition.
                                 0: with linear elements
                                 1: with quadratic elemens.
        */
        NodePartitionedMesh(const std::string &name,
                            const std::vector<Node*> &nodes,
                            const std::vector<unsigned> &glb_node_ids,
                            const std::vector<Element*> &elements,
                            const std::vector<short*> &g_elem_data,
                            const std::size_t n_nghost_elem,
                            const unsigned nnodes_global[],
                            const unsigned nnodes_active[])
            : Mesh(name, nodes, elements, false),
              _global_node_ids(glb_node_ids),
              _n_nghost_elem(n_nghost_elem),
              _nnodes_global {nnodes_global[0], nnodes_global[1] },
                         _nnodes_active {nnodes_active[0], nnodes_active[1] },
                         _act_nodes_ids_of_ghost_element(g_elem_data)
        {
        }

        ~NodePartitionedMesh()
        {
            for(auto ele_data : _act_nodes_ids_of_ghost_element)
            {
                delete [] ele_data;
            }
        }

        /// Get the number of nodes of the global mesh for linear elements.
        unsigned getNGlobalBaseNodes() const
        {
            return _nnodes_global[0];
        }

        /// Get the number of all nodes of the global mesh.
        unsigned getNGlobalNodes() const
        {
            return _nnodes_global[1];
        }

        /// Get the number of the active nodes of the partition for linear elements.
        unsigned getNActiveBaseNodes() const
        {
            return _nnodes_active[0];
        }

        /// Get the number of all active nodes of the partition.
        unsigned getNActiveNodes() const
        {
            return _nnodes_active[1];
        }

        /*!
          \brief Get the number of active nodes of a ghost element for linear interpolation
          \param gelem_id Index of ghost element
          \return         The first or the second element of array by
                           _nnodes_active[static_cast<gelem_id> , which stores
                          the number of active nodes either for linear or high
                          order element of an ghost element.
        */
        short getNGhostElementActiveBaseNodes(const unsigned gelem_id) const
        {
            return _act_nodes_ids_of_ghost_element[gelem_id][0];
        }

        /*!
          \brief Get the number of all active nodes of a ghost element
          \param gelem_id Index of ghost element
          \return         The first or the second element of array by
                           _nnodes_active[static_cast<gelem_id> , which stores
                          the number of active nodes either for linear or high
                          order element of an ghost element.
        */
        short getNGhostElementActiveNodes(const unsigned gelem_id) const
        {
            return _act_nodes_ids_of_ghost_element[gelem_id][1];
        }


        /*!
         \brief Get local IDs of the active nodes of a ghost element by a pointer to
                an array.
         \param gelem_id Index of ghost element.
        */
        short *getGhostElementActiveNodes(const unsigned gelem_id) const
        {
            return &_act_nodes_ids_of_ghost_element[gelem_id][2];
        }

        /// Get the largest ID of active nodes for higher order elements in a partition.
        unsigned getLargestActiveNodeID() const
        {
            // Note: _nodes.size() should be changed once the high order element is condidered
            // in the root class.
            return static_cast<unsigned>( _nodes.size() ) + _nnodes_active[1] - _nnodes_active[0];
        }

        /// Get the number of non-ghost elements, or the start entry ID of ghost elements in element vector.
        size_t getNNonGhostElements() const
        {
            return _n_nghost_elem;
        }

    private:
        /// Global IDs of nodes of a partition
        std::vector<unsigned> _global_node_ids;

        /// Number of non-ghost elements, or the ID of the start entry of ghost elements in _elements vector.
        std::size_t _n_nghost_elem;

        /// Number of nodes of the whole mesh. 0: for linear elements; 1: for quadratic elements.
        unsigned _nnodes_global[2];

        /// Number of the active nodes. 0: for linear elements; 1: for quadratic elements.
        unsigned  _nnodes_active[2];

        /*! Active node indices of each ghost elements.
            In each element of the vector, an integer array, the first and the second element of the array
            stores the numbers of active nodes either for linear and high order element, respectively,
            while the remaining elements of the array are for local IDs of active nodes.
        */
        std::vector<short *> _act_nodes_ids_of_ghost_element;

        friend FileIO::readNodePartitionedMesh;
};

} // end of namespace

#endif // end of #ifndef NODE_PARTITIONED_MESH_H_

