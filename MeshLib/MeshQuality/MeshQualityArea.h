/*
 * MeshQualityArea.h
 *
 * 2011/03/17 KR Initial Implementation
 */

#ifndef MESHQUALITYAREA_H_
#define MESHQUALITYAREA_H_

#include "MeshQualityChecker.h"

namespace MeshLib
{
class MeshQualityArea : public MeshQualityChecker
{
public:
	MeshQualityArea(Mesh const* const mesh);
	virtual ~MeshQualityArea() {}

	virtual void check ();
};
}

#endif /* MESHQUALITYAREA_H_ */
