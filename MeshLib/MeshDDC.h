/*!
  \file MeshDDC.h
  \author Wenqing Wang
  \date   2014.06
  \brief  Definition of the Mesh class for a partition for parallel computing within the
          framework of domain decomposition (DDC).

  \copyright
  Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

 */

#ifndef MESH_DDC_H_
#define MESH_DDC_H_

#include <cstdlib>
#include <string>
#include <vector>

#include "Mesh"

namespace MeshLib
{
class Node;
class Element;

/// A subdomain mesh.
class MeshDDC : public Mesh
{
    public:
        /// Constructor using a mesh name and an array of nodes and elements
        MeshDDC(const std::string &name,
                const std::vector<Node*> &nodes,
                const std::vector<Element*> &elements,
                const std::vector<unsigned> g_node_ids,
                const std::vector<unsigned*> l_act_nodes_ids);

        /*!
            Get global index of a node
            \param sdom_node_id Subdmain node ID
        */
        const unsigned getGlobalNodeID(const unsigned sdom_node_id);

        /*!
            Get number of active nodes of an element
            \param elem_id Element ID
            \param order   Order of element, could be 0 or 1
        */
        const unsigned getNumberActiveNodes(const unsigned elem_id, const int order)
        {
            _l_act_nodes_ids[elem_id][order];
        }

        /*!
            Get local indicies of active nodes
        */
        const unsigned *getLocalActiveNodeIDs(const unsigned elem_id);
        {
            return &_l_act_nodes_ids[elem_id][2];
        }

    private:
        /// Global indicies of nodes. the vector size is equal to the number of all nodes
        std::vector<unsigned> _g_node_ids;

        /// Local indices of active nodes of each elements.
        /// The vector size is equal to the number of elements
        std::vector<unsigned*> _l_act_nodes_ids;

};

} // end of namespace

#endif // end of #ifndef MESH_DDC_H_

