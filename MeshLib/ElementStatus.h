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
 * Manages active/inactive mesh elements and their nodes
 */
class ElementStatus
{

public:
	/// Constructor
	ElementStatus(MeshLib::Mesh const*const mesh);

	/// Returns a vector of active element IDs
	std::vector<unsigned> getActiveElements() const;

	/// Returns a vector of active node IDs
	std::vector<unsigned> getActiveNodes() const;
	
	/// Returns the status of element i
	bool getElementStatus(unsigned i) const { return _element_status[i]; }

	/// Returns a vector of active elements connected to a node
	std::vector<MeshLib::Element*> getActiveElementsAtNode(unsigned node_id) const;
	
	/// Returns the total number of active nodes
	unsigned getNActiveNodes() const;

	/// Returns the total number of active elements
	unsigned getNActiveElements() const;

	/// Returns the status of element 
	bool isActive(unsigned i) const { return _element_status[i]; } 

	/// Sets the status of element i
	void setElementStatus(unsigned i, bool status);

	/// Sets the status of material group i
	void setMaterialStatus(unsigned material_id, bool status);

	~ElementStatus() {};

protected:
	MeshLib::Mesh const*const _mesh;
	std::vector<bool> _element_status;
	std::vector<char> _active_nodes;

}; /* class */

} /* namespace */

#endif /* ELEMENTSTATUS_H_ */

