/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file FemNode.h
 *
 * Created on 2012-05-02 by Karsten Rink
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
	FemNode(double const*const coords, unsigned id = std::numeric_limits<unsigned>::max());

	/// Constructor using single coordinates
	FemNode(double x, double y, double z, unsigned id = std::numeric_limits<unsigned>::max());

	/// Constructor using a mesh node
	FemNode(const Node &node);

	/// Copy constructor
	FemNode(const FemNode &node);

	/// Destructor
	virtual ~FemNode();

}; /* class */

} /* namespace */

#endif /* FEMNODE_H_ */

