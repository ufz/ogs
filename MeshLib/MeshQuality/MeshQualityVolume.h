/**
 * \file
 * \author Thomas Fischer
 * \date   2011-03-03
 * \brief  Definition of the MeshQualityVolume class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
