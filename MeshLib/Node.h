/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Node class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NODE_H_
#define NODE_H_

#include <cstdlib>
#include <limits>
#include <vector>
#include <set>

#include "PointWithID.h"
#include "Mesh.h"
#include "MeshEditing/removeMeshNodes.h"
#include "MeshSurfaceExtraction.h"
#include "MeshGenerators/MeshLayerMapper.h"

namespace MeshLib {

class Element;

/**
 * A mesh node with coordinates in 3D space.
 */
class Node : public GeoLib::PointWithID
{
	/* friend functions: */
	friend bool MeshLayerMapper::LayerMapping(MeshLib::Mesh &mesh, const GeoLib::Raster &raster, const unsigned nLayers,
		                                      const unsigned layer_id, double noDataReplacementValue);
	friend MeshLib::Mesh* MeshLayerMapper::blendLayersWithSurface(MeshLib::Mesh &mesh, const unsigned nLayers, const std::string &dem_raster);
	friend MeshLib::Mesh* MeshSurfaceExtraction::getMeshSurface(const MeshLib::Mesh &mesh, const MathLib::Vector3 &dir, bool keep3dMeshIds);

	/* friend classes: */
	friend class Mesh;
	friend class MeshRevision;


public:
	/// Constructor using a coordinate array
	Node(const double coords[3], unsigned id = std::numeric_limits<unsigned>::max());

	/// Constructor using single coordinates
	Node(double x, double y, double z, unsigned id = std::numeric_limits<unsigned>::max());

	/// Copy constructor
	Node(const Node &node);

	/// Return all the nodes connected to this one
	const std::vector<MeshLib::Node*>& getConnectedNodes() const { return _connected_nodes; }

	/// Get an element the node is part of.
	const Element* getElement(unsigned idx) const { return _elements[idx]; }

	/// Get all elements the node is part of.
	const std::vector<Element*>& getElements() const { return _elements; }

	/// Get number of elements the node is part of.
	std::size_t getNElements() const { return _elements.size(); }

	/// Destructor
	virtual ~Node();

protected:
	/**
	 * Add an element the node is part of.
	 * This method is called by Mesh::addElement(Element*), see friend definition.
	 */
	void addElement(Element* elem) { _elements.push_back(elem); }

	void setConnectedNodes(std::vector<Node*> &connected_nodes) { this->_connected_nodes = connected_nodes; }

	/// Sets the ID of a node to the given value.
	void setID(unsigned id) { this->_id = id; }

	/// Update coordinates of a node.
	/// This method automatically also updates the areas/volumes of all connected elements.
	virtual void updateCoordinates(double x, double y, double z);

	std::vector<Node*> _connected_nodes;
	std::vector<Element*> _elements;

}; /* class */

} /* namespace */

#endif /* NODE_H_ */

