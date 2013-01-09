/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file GMSHAdaptiveMeshDensity.h
 *
 *  Created on 2012-03-05 by Thomas Fischer
 */

#ifndef GMSHADAPTIVEMESHDENSITY_H_
#define GMSHADAPTIVEMESHDENSITY_H_

// FileIO
#include "GMSHMeshDensityStrategy.h"

// GeoLib
#include "Point.h"
#include "QuadTree.h"

namespace GeoLib
{
class Polygon;
#ifndef NDEBUG
class Polyline;
#endif
}

namespace FileIO
{
class GMSHAdaptiveMeshDensity : public GMSHMeshDensityStrategy
{
public:
	GMSHAdaptiveMeshDensity(double pnt_density,
	                        double station_density,
	                        std::size_t max_pnts_per_leaf);
	virtual ~GMSHAdaptiveMeshDensity();
	void init(std::vector<GeoLib::Point const*> const& pnts);
	double getMeshDensityAtPoint(GeoLib::Point const* const pnt) const;
	void addPoints(std::vector<GeoLib::Point const*> const& pnts);
	double getMeshDensityAtStation(GeoLib::Point const* const) const;
	void getSteinerPoints (std::vector<GeoLib::Point*> & pnts,
	                       std::size_t additional_levels = 0) const;
#ifndef NDEBUG
	void getQuadTreeGeometry(std::vector<GeoLib::Point*> &pnts,
	                         std::vector<GeoLib::Polyline*> &plys) const;
#endif

private:
	double _pnt_density;
	double _station_density;
	std::size_t _max_pnts_per_leaf;
	GeoLib::QuadTree<GeoLib::Point> *_quad_tree;
};
} // end namespace FileIO

#endif /* GMSHADAPTIVEMESHDENSITY_H_ */
