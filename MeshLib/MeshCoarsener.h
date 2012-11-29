/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file MeshCoarsener.h
 *
 *  Created on  Aug 3, 2012 by Thomas Fischer
 */

#ifndef MESHCOARSENER_H_
#define MESHCOARSENER_H_

#include <vector>
#include <cstddef>

// forward declaration
namespace MeshLib {
class Mesh;
}

namespace MeshLib {

/**
 * The class MeshCoarsener merges mesh nodes that have a smaller
 * distance than a (user) given minimal distance.
 */
class MeshCoarsener {
public:
	/**
	 * Constructor of class MeshCoarsener that takes the mesh object that should be coarsened.
	 * @param mesh the mesh object
	 */
	MeshCoarsener(Mesh const*const mesh);

	/**
	 * destructor.
	 */
	virtual ~MeshCoarsener();

	/**
	 * create new mesh and apply the coarsening process to the mesh
	 * @param min_distance
	 */
	Mesh* operator() (double min_distance);

private:
	Mesh const*const _orig_mesh;
};

}

#endif /* MESHCOARSENER_H_ */
