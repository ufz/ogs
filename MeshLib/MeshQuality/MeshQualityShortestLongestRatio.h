/*
 * MeshQualityShortestLongestRatio.h
 *
 *  Created on: Mar 3, 2011
 *      Author: TF
 */

#ifndef MESHQUALITYSHORTESTLONGESTRATIO_H_
#define MESHQUALITYSHORTESTLONGESTRATIO_H_

#include "MeshQualityChecker.h"
#include "Point.h"

namespace MeshLib
{
class MeshQualityShortestLongestRatio : public MeshQualityChecker
{
public:
	MeshQualityShortestLongestRatio(Mesh const* const mesh);
	virtual ~MeshQualityShortestLongestRatio () {}

	virtual void check ();

private:
	double checkTriangle (GeoLib::Point const* const a,
	                      GeoLib::Point const* const b,
	                      GeoLib::Point const* const c) const;
	double checkQuad (GeoLib::Point const* const a,
	                  GeoLib::Point const* const b,
	                  GeoLib::Point const* const c,
	                  GeoLib::Point const* const d) const;
	double checkTetrahedron (GeoLib::Point const* const a,
	                         GeoLib::Point const* const b,
	                         GeoLib::Point const* const c,
	                         GeoLib::Point const* const d) const;
	double checkPrism (std::vector<const GeoLib::Point*> const & pnts) const;
	double checkHexahedron (std::vector<const GeoLib::Point*> const & pnts) const;
};
}

#endif /* MESHQUALITYSHORTESTLONGESTRATIO_H_ */
