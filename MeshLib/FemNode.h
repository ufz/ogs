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

class FemNode : public Node
{
public:
	FemNode(double const*const coords, size_t id = std::numeric_limits<size_t>::max());
	FemNode(double x, double y, double z, size_t id = std::numeric_limits<size_t>::max());
	FemNode(const Node &node);
	FemNode(const FemNode &node);
	~FemNode();

}; /* class */

} /* namespace */

#endif /* FEMNODE_H_ */
