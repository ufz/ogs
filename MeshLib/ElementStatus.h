/**
 * \file   ElementStatus.h
 * \author Karsten Rink
 * \date   2012-12-18
 * \brief  Definition of the ElementStatus class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENTSTATUS_H_
#define ELEMENTSTATUS_H_

#include <vector>

namespace MeshLib {
    class Mesh;
    class Element;
    class Node;

class ElementStatus
{

public:
    /// Constructor assuming all nodes and elements
    explicit ElementStatus(MeshLib::Mesh const* const mesh,
                           bool hasAnyInactive = false);

    /// Constructor taking a vector of inactive material IDs
    ElementStatus(MeshLib::Mesh const* const mesh,
                  std::vector<int> const& vec_inactive_matIDs);

    /// Returns a vector of active element IDs
    std::vector<MeshLib::Element*> const& getActiveElements() const;

    /// Returns a vector of active node IDs
    std::vector<MeshLib::Node*> const& getActiveNodes() const;

    /// Returns the status of element i
    bool isActive(std::size_t i) const { return _element_status[i]; }

    /// Returns the status of the given node
    bool isActiveNode(MeshLib::Node const* node) const;

    /// Returns a vector of active elements connected to a node
    std::vector<MeshLib::Element*> getActiveElementsAtNode(std::size_t node_id) const;

    /// Returns the total number of active nodes
    std::size_t getNumberOfActiveNodes() const;

    /// Returns the total number of active elements
    std::size_t getNumberOfActiveElements() const;

protected:
    /// Sets the status of element i
    void setElementStatus(std::size_t i, bool status);

    /// The mesh for which the element status is administrated
    MeshLib::Mesh const*const _mesh;
    /// Element status for each mesh element (active/inactive = true/false)
    std::vector<bool> _element_status;
    /// Node status for each mesh node (value = number of active elements connected to node, 0 means inactive)
    std::vector<unsigned char> _active_nodes;

    bool const _hasAnyInactive;
    std::vector<MeshLib::Node*> _vec_active_nodes;
    std::vector<MeshLib::Element*> _vec_active_eles;

}; /* class */

} /* namespace */

#endif /* ELEMENTSTATUS_H_ */
