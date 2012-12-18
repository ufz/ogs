/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file ElementStatus.h
 *
 * Created on 2012-12-18 by Karsten Rink
 */

#ifndef ELEMENTSTATUS_H_
#define ELEMENTSTATUS_H_

#include "Mesh.h"

namespace MeshLib {
	
/**
 * An object to store which elements of the mesh are active and which are inactive
 */
class ElementStatus
{

public:
	/// Constructor using mesh
	ElementStatus(Mesh const*const mesh);

	bool getElementStatus(unsigned i) { return _status[i]; };
	
	bool isActive(unsigned i) { return _status[i]; };

	void setActive(unsigned i) { _status[i] = true; };

	void setInactive(unsigned i) { _status[i] = false; };

	void setElementStatus(unsigned i, bool status) { _status[i] = status; };

	~ElementStatus() {};

protected:
	Mesh const*const _mesh;
	std::vector<bool> _status;

}; /* class */

} /* namespace */

#endif /* ELEMENTSTATUS_H_ */

