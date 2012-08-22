/*
 *  Created on: Mar 3, 2011
 *      Author: TF
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
