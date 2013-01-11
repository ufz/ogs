/**
 * \file
 * \author Thomas Fischer
 * \date   2011-03-03
 * \brief  Definition of the MeshQualityShortestLongestRatio class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
