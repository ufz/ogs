/**
 * Mesh.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef MESH_H_
#define MESH_H_

#include <cstdlib>
#include <string>
#include <vector>

namespace MeshLib {

class Node;
class Element;

/**
 * A basic mesh.
 */
class Mesh
{
public:
	/// Constructor using a mesh name and an array of nodes and elements
	Mesh(const std::string &name, const std::vector<Node*> &nodes, const std::vector<Element*> &elements);

	/// Copy constructor
	Mesh(const Mesh &mesh);

	/// Constructor for reading a mesh from a file
	Mesh(const std::string &file_name);

	/// Destructor
	virtual ~Mesh();

	/// Add a node to the mesh.
	void addNode(Node* node);

	/// Add an element to the mesh.
	void addElement(Element* elem);

	/// Get the node with the given index.
	const Node* getNode(size_t idx) const { return _nodes[idx]; };

	/// Get the element with the given index.
	const Element* getElement(size_t idx) const { return _elements[idx]; };

	/// Get the number of elements
	size_t getNElements() const { return _elements.size(); };

	/// Get the number of nodes
	size_t getNNodes() const { return _nodes.size(); };

	/// Get name of the mesh.
	const std::string getName() const { return _name; };

	/// Get the nodes-vector for the mesh.
	const std::vector<Node*> getNodes() const { return _nodes; };

	/// Get the element-vector for the mesh.
	const std::vector<Element*> getElements() const { return _elements; };

protected:
	/// Checks the coordinates of all mesh nodes and removes identical nodes. Elements are adapted accordingly.
	void removeIdenticalNodes();

	std::string _name;
	std::vector<Node*> _nodes;
	std::vector<Element*> _elements;

}; /* class */

} /* namespace */

#endif /* MESH_H_ */

