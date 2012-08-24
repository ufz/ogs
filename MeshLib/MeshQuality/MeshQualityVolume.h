/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file MeshQualityVolume.h
 *
 *  Created on 2011-03-03 by Thomas Fischer
 */

#ifndef MESHQUALITYVOLUME_H_
#define MESHQUALITYVOLUME_H_

#include "MeshQualityChecker.h"

namespace MeshLib
{
class MeshQualityVolume : public MeshQualityChecker
{
public:
	MeshQualityVolume(Mesh const* const mesh);
	virtual ~MeshQualityVolume() {}

	virtual void check ();
};
}

#endif /* MESHQUALITYVOLUME_H_ */
