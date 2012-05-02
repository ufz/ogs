/**
 * FemNode.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef FEMNODE_H_
#define FEMNODE_H_

#include "Node.h"

namespace MeshLib {

/**
 * A mesh node for finite element meshes.
 */
class FemNode : public Node
{
public:
	/// Constructor using a coordinate array
	FemNode(double const*const coords, size_t id = std::numeric_limits<size_t>::max());
	
	/// Constructor using single coordinates
	FemNode(double x, double y, double z, size_t id = std::numeric_limits<size_t>::max());

	/// Constructor using a mesh node
	FemNode(const Node &node);

	/// Copy constructor
	FemNode(const FemNode &node);

	/// Destructor
	~FemNode();

}; /* class */

} /* namespace */

#endif /* FEMNODE_H_ */
