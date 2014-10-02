/*!
  \file NodePartitionedMesh.h
  \author Wenqing Wang
  \date   2014.06
  \brief  Definition of mesh class for partitioned mesh (by node) for parallel computing within the
          framework of domain decomposition (DDC).

  \copyright
  Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
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

/// Element order for reading partitioned mesh data
enum class ElementOrder
{
    LINEAR = 0,
    QUADRATIC = 1,
};

/// A subdomain mesh.
class NodePartitionedMesh : public Mesh
{
    public:
        /*!
            \brief  Constructor
            \param name          Name assigned to the mesh.
            \param nodes         Vector for nodes.
            \param elements      Vector for elements. Ghost elements are stored
                                 after regular elements.
            \param g_elem_data   Element wise local node IDs of active nodes.
            \param start_id_gele Start ID of the entry of ghost element in
                                 the element vector.
            \param nnodes_global Number of nodes of the global mesh.
                                 0: with linear elements
                                 1: with quadratic elemens.
            \param nnodes_active Number of active nodes of the partition.
                                 0: with linear elements
                                 1: with quadratic elemens.
        */
        NodePartitionedMesh(const std::string &name,
                            const std::vector<Node*> &nodes,
                            const std::vector<Element*> &elements,
                            const std::vector<short*> &g_elem_data,
                            const std::size_t start_id_g_elem,
                            const unsigned nnodes_global[],
                            const unsigned nnodes_active[])
            : Mesh(name, nodes, elements, false),
              _start_id_g_elem(start_id_g_elem),
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

        /*!
            \brief Get the number of nodes of the whole mesh.
            \param order The order of elements (0 or 1).
                         Its default value is 0 for linear elements.
        */
        unsigned getNGlobalNodes(const ElementOrder order = ElementOrder::LINEAR) const
        {
            return _nnodes_global[static_cast<int>(order)];
        }

        /*!
            \brief Get the number of the active nodes of the partition.
            \param order The order of elements (0 or 1).
                         Its default value is 0 for linear elements.
        */
        unsigned getNActiveNodes(const ElementOrder order = ElementOrder::LINEAR) const
        {
            return _nnodes_active[static_cast<int>(order)];
        }

        /*!
          \brief Get the number of active nodes of a ghost element
          \param gelem_id Index of ghost element
          \param order    The order of elements (either 0 or 1).
          \return         The first or the second element of array by
                           _nnodes_active[static_cast<gelem_id> , which stores
                          the number of active nodes either for linear or high
                          order element of an ghost element.
        */
        unsigned getNGhostElementActiveNodes(const unsigned gelem_id,
                                             const ElementOrder order = ElementOrder::LINEAR) const
        {
            return _act_nodes_ids_of_ghost_element[gelem_id][static_cast<int>(order)];
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

        /// Get the start entry ID of ghost elements in element vector.
        size_t getStartIndexOfGhostElement() const
        {
            return _start_id_g_elem;
        }

    private:
        /// ID of the start entry of ghost elements in _elements vector.
        std::size_t _start_id_g_elem;

        /// Number of nodes of the whole mesh. 0: for linear elements; 1: for quadratic elements.
        unsigned _nnodes_global[2];

        /// Number of the active nodes. 0: for linear elements; 1: for quadratic elements.
        unsigned  _nnodes_active[2];

        /*! Active node indices of each ghost elements.
            In each element of the vector, an integger array, the first and the second element of the array
            stores the numbers of active nodes either for linear and high order element, repectively,
            while the remaining elements of the array are for local IDs of active .
        */
        std::vector<short *> _act_nodes_ids_of_ghost_element;

        friend FileIO::readNodePartitionedMesh;
};

} // end of namespace

#endif // end of #ifndef NODE_PARTITIONED_MESH_H_

