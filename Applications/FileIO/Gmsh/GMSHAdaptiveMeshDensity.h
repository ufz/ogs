/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GMSHADAPTIVEMESHDENSITY_H_
#define GMSHADAPTIVEMESHDENSITY_H_

#include "GMSHMeshDensityStrategy.h"

namespace GeoLib
{
class Polygon;
class Point;
template <typename T> class QuadTree;
#ifndef NDEBUG
class Polyline;
#endif
}

namespace FileIO
{
namespace GMSH
{

class GMSHAdaptiveMeshDensity : public GMSHMeshDensityStrategy
{
public:
    GMSHAdaptiveMeshDensity(double pnt_density,
                            double station_density,
                            std::size_t max_pnts_per_leaf);
    virtual ~GMSHAdaptiveMeshDensity();
    void initialize(std::vector<GeoLib::Point const*> const& pnts);
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

}
} // end namespace FileIO

#endif /* GMSHADAPTIVEMESHDENSITY_H_ */
