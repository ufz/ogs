/*
 * GMSHAdaptiveMeshDensity.h
 *
 *  Created on: Mar 5, 2012
 *      Author: TF
 */

#ifndef GMSHADAPTIVEMESHDENSITY_H_
#define GMSHADAPTIVEMESHDENSITY_H_

// FileIO
#include "GMSHMeshDensityStrategy.h"

// GeoLib
#include "Point.h"
#include "QuadTree.h"
#ifndef NDEBUG
#include "Polyline.h"
#endif

namespace GeoLib {
class Polygon;
}

namespace FileIO {

class GMSHAdaptiveMeshDensity: public GMSHMeshDensityStrategy {
public:
	GMSHAdaptiveMeshDensity(double pnt_density, double station_density, size_t max_pnts_per_leaf);
	virtual ~GMSHAdaptiveMeshDensity();
	void init(std::vector<GeoLib::Point const*> const& pnts);
	double getMeshDensityAtPoint(GeoLib::Point const*const pnt) const;
	void addPoints(std::vector<GeoLib::Point const*> const& pnts);
	double getMeshDensityAtStation(GeoLib::Point const*const) const;
	void getSteinerPoints (std::vector<GeoLib::Point*> & pnts, size_t additional_levels = 0) const;
#ifndef NDEBUG
	void getQuadTreeGeometry(std::vector<GeoLib::Point*> &pnts, std::vector<GeoLib::Polyline*> &plys) const;
#endif

private:
	double _pnt_density;
	double _station_density;
	size_t _max_pnts_per_leaf;
	GeoLib::QuadTree<GeoLib::Point> *_quad_tree;
};

} // end namespace FileIO

#endif /* GMSHADAPTIVEMESHDENSITY_H_ */
