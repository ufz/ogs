/**
 * \file ElementDDC
 * \author Wenqing Wang
 * \date   2014-06
 * \brief  Definition of the Element class for domain decomposition (DDC) approach.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENT_DDC_H_
#define ELEMENT_DDC_H_

#include "Element.h"

namespace MeshLib
{


/**
 * Virtual base class for mesh elements in subdomain.
 */
class ElementDDC : puclic Element
{
    public:
        ElementDDC(const unsigned *active_node_ids);

        /*!
              Get number of active nodes
             \param order   Order of element, could be 0 or 1
         */
        virtual unsigned getNNodes(bool all = false) const
        {
            _act_nodes_num[all];
        }

        /*!
            Get local indicies of active nodes
            \param local_id ID of entry of node pointer array
        */
        virtual const unsigned getLocalActiveNodeID(const unsigned local_id);
        {
            return _l_act_nodes_ids[local_id];
        }

    private:
        /*!
           Number of active nodes
           [0]: number of active nodes
           [1]: number of active nodes for a high order element
        */
        unsigned _act_nodes_num;
        /*!
           Local indices of active nodes of each elements.
        */
        unsigned* _l_act_nodes_ids;

} // end of namespace

#endif // end of #ifndef ELEMENT_DDC_H_ 

