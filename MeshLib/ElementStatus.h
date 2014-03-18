/**
 * \file   ElementStatus.h
 * \author Karsten Rink
 * \date   2012-12-18
 * \brief  Definition of the ElementStatus class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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

/**
 * Manages the active/inactive flags for mesh elements and their nodes
 */
class ElementStatus
{

public:
	/// Constructor
	explicit ElementStatus(MeshLib::Mesh const*const mesh);

	/// Returns a vector of active element IDs
	std::vector<std::size_t> getActiveElements() const;

	/// Returns a vector of active node IDs
	std::vector<std::size_t> getActiveNodes() const;
	
	/// Returns the status of element i
	bool getElementStatus(std::size_t i) const { return _element_status[i]; }

	/// Returns a vector of active elements connected to a node
	std::vector<std::size_t> getActiveElementsAtNode(std::size_t node_id) const;
	
	/// Returns the total number of active nodes
	std::size_t getNActiveNodes() const;

	/// Returns the total number of active elements
	std::size_t getNActiveElements() const;

	/// Activates/Deactives all mesh elements
	void setAll(bool status);

	/// Sets the status of element i
	void setElementStatus(std::size_t i, bool status);

	/// Sets the status of material group i
	void setMaterialStatus(unsigned material_id, bool status);

	~ElementStatus() {};

protected:
	/// The mesh for which the element status is administrated
	MeshLib::Mesh const*const _mesh;
	/// Element status for each mesh element (active/inactive = true/false)
	std::vector<bool> _element_status;
	/// Node status for each mesh node (value = number of active elements connected to node, 0 means inactive)
	std::vector<unsigned char> _active_nodes;

}; /* class */

} /* namespace */

#endif /* ELEMENTSTATUS_H_ */

