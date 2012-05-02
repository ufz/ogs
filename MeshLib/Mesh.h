/**
 * Mesh.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef MESH_H_
#define MESH_H_

#include <cstdlib>
#include <vector>

namespace MeshLib {

class Node;
class Element;

class Mesh
{
public:
	Mesh(const std::string &name, const std::vector<Node*> &nodes, const std::vector<Element*> &elements);
	Mesh(const Mesh &mesh);
	Mesh(const std::string &file_name);
	virtual ~Mesh();

	const std::string getName() const { return _name; };
	const std::vector<Node*> getNodes() const { return _nodes; };
	const std::vector<Element*> getElements() const { return _elements; };

protected:
	void removeIdenticalNodes();

	std::string _name;
	std::vector<Node*> _nodes;
	std::vector<Element*> _elements;

}; /* class */

} /* namespace */

#endif /* MESH_H_ */
