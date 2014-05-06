/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <fstream>
#include <memory>
#include <vector>

#include <logog/include/logog.hpp>

#include "BaseLib/BuildInfo.h"
#include "BaseLib/FileTools.h"

#include "Applications/FileIO/Gmsh/GMSHInterface.h"
#include "Applications/FileIO/Gmsh/GMSHAdaptiveMeshDensity.h"
#include "Applications/FileIO/Gmsh/GMSHFixedMeshDensity.h"
#include "Applications/FileIO/Gmsh/GMSHPoint.h"
#include "Applications/FileIO/Gmsh/GMSHPolygonTree.h"
#include "Applications/FileIO/Gmsh/GMSHMeshDensityStrategy.h"

#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polygon.h"
#include "GeoLib/PolylineWithSegmentMarker.h"
#include "GeoLib/PolygonWithSegmentMarker.h"

namespace FileIO
{
namespace GMSH
{
GMSHInterface::GMSHInterface(GeoLib::GEOObjects& geo_objs,
                             bool /*include_stations_as_constraints*/,
                             GMSH::MeshDensityAlgorithm mesh_density_algorithm,
                             double param1, double param2, std::size_t param3,
                             std::vector<std::string>& selected_geometries,
                             bool rotate, bool keep_preprocessed_geometry)
    : _n_lines(0),
      _n_plane_sfc(0),
      _geo_objs(geo_objs),
      _selected_geometries(selected_geometries),
      _rotate(rotate),
      _keep_preprocessed_geometry(keep_preprocessed_geometry)
{
    switch (mesh_density_algorithm) {
    case GMSH::MeshDensityAlgorithm::FixedMeshDensity:
        _mesh_density_strategy = new GMSH::GMSHFixedMeshDensity(param1);
        break;
    case GMSH::MeshDensityAlgorithm::AdaptiveMeshDensity:
        _mesh_density_strategy = new GMSH::GMSHAdaptiveMeshDensity(param1, param2, param3);
        break;
    }
}

GMSHInterface::~GMSHInterface()
{
    for (auto * gmsh_pnt : _gmsh_pnts)
        delete gmsh_pnt;
    delete _mesh_density_strategy;
    for (auto * polygon_tree : _polygon_tree_list)
        delete polygon_tree;
}

bool GMSHInterface::write()
{
    _out << "// GMSH input file created by OpenGeoSys " << BaseLib::BuildInfo::git_describe;
#ifdef BUILD_TIMESTAMP
    _out << " built on " << BaseLib::BuildInfo::build_timestamp;
#endif
    _out << "\n\n";

    return writeGMSHInputFile(_out) <= 0;
}

int GMSHInterface::writeGMSHInputFile(std::ostream& out)
{
    DBUG("GMSHInterface::writeGMSHInputFile(): get data from GEOObjects.");

    if (_selected_geometries.empty())
        return 1;

    // *** get and merge data from _geo_objs
    if (_selected_geometries.size() > 1) {
        _gmsh_geo_name = "GMSHGeometry";
        if (_geo_objs.mergeGeometries(_selected_geometries, _gmsh_geo_name))
            return 2;
    } else {
        _gmsh_geo_name = _selected_geometries[0];
        _keep_preprocessed_geometry = true;
    }

    auto* merged_pnts(const_cast<std::vector<GeoLib::Point*>*>(
        _geo_objs.getPointVec(_gmsh_geo_name)));
    if (! merged_pnts) {
        ERR("GMSHInterface::writeGMSHInputFile(): Did not found any points.");
        return 2;
    }

    if (_rotate) {
        // Rotate points to the x-y-plane.
        _inverse_rot_mat = GeoLib::rotatePointsToXY(*merged_pnts);
        // Compute inverse rotation matrix to reverse the rotation later on.
        _inverse_rot_mat.transposeInPlace();
    } else {
        // project data on the x-y-plane
        _inverse_rot_mat(0,0) = 1.0;
        _inverse_rot_mat(1,1) = 1.0;
        _inverse_rot_mat(2,2) = 1.0;
        for (auto pnt : *merged_pnts)
            (*pnt)[2] = 0.0;
    }

    std::vector<GeoLib::Polyline*> const* merged_plys(
        _geo_objs.getPolylineVec(_gmsh_geo_name));
    DBUG("GMSHInterface::writeGMSHInputFile(): Obtained data.");

    if (!merged_plys) {
        ERR("GMSHInterface::writeGMSHInputFile(): Did not find any polylines.");
        return 2;
    }

    // *** compute and insert all intersection points between polylines
    GeoLib::PointVec& pnt_vec(*const_cast<GeoLib::PointVec*>(
        _geo_objs.getPointVecObj(_gmsh_geo_name)));
    GeoLib::computeAndInsertAllIntersectionPoints(pnt_vec,
        *(const_cast<std::vector<GeoLib::Polyline*>*>(merged_plys)));

    // *** compute topological hierarchy of polygons
    for (auto polyline : *merged_plys) {
        if (!polyline->isClosed()) {
            continue;
        }
        _polygon_tree_list.push_back(new GMSH::GMSHPolygonTree(
            new GeoLib::PolygonWithSegmentMarker(*polyline), nullptr, _geo_objs,
            _gmsh_geo_name, _mesh_density_strategy));
    }
    DBUG(
        "GMSHInterface::writeGMSHInputFile(): Computed topological hierarchy - "
        "detected %d polygons.",
        _polygon_tree_list.size());
    GeoLib::createPolygonTrees<GMSH::GMSHPolygonTree>(_polygon_tree_list);
    DBUG(
        "GMSHInterface::writeGMSHInputFile(): Computed topological hierarchy - "
        "calculated %d polygon trees.",
        _polygon_tree_list.size());

    // *** Mark in each polygon tree the segments shared by two polygons.
    for (auto* polygon_tree : _polygon_tree_list)
    {
        polygon_tree->markSharedSegments();
    }

    // *** insert stations and polylines (except polygons) in the appropriate object of
    //     class GMSHPolygonTree
    // *** insert stations
    auto gmsh_stations = std::make_unique<std::vector<GeoLib::Point*>>();
    for (auto const& geometry_name : _selected_geometries) {
        auto const* stations(_geo_objs.getStationVec(geometry_name));
        if (stations) {
            for (auto * station : *stations) {
                bool found(false);
                for (auto it(_polygon_tree_list.begin());
                    it != _polygon_tree_list.end() && !found; ++it) {
                    gmsh_stations->emplace_back(new GeoLib::Station(
                        *static_cast<GeoLib::Station*>(station)));
                    if ((*it)->insertStation(gmsh_stations->back())) {
                        found = true;
                    }
                }
            }
        }
    }
    std::string gmsh_stations_name(_gmsh_geo_name+"-Stations");
    if (! gmsh_stations->empty()) {
        _geo_objs.addStationVec(std::move(gmsh_stations), gmsh_stations_name);
    }

    // *** insert polylines
    for (auto polyline : *merged_plys) {
        if (!polyline->isClosed()) {
            for (auto * polygon_tree : _polygon_tree_list) {
                auto polyline_with_segment_marker =
                    new GeoLib::PolylineWithSegmentMarker(*polyline);
                polygon_tree->insertPolyline(polyline_with_segment_marker);
            }
        }
    }

    // *** init mesh density strategies
    for (auto& polygon_tree : _polygon_tree_list)
    {
        polygon_tree->initMeshDensityStrategy();
    }

    // *** create GMSH data structures
    const std::size_t n_merged_pnts(merged_pnts->size());
    _gmsh_pnts.resize(n_merged_pnts);
    for (std::size_t k(0); k<n_merged_pnts; k++) {
        _gmsh_pnts[k] = nullptr;
    }
    for (auto& polygon_tree : _polygon_tree_list)
    {
        polygon_tree->createGMSHPoints(_gmsh_pnts);
    }

    // *** finally write data :-)
    writePoints(out);
    std::size_t pnt_id_offset(_gmsh_pnts.size());
    for (auto* polygon_tree : _polygon_tree_list)
    {
        polygon_tree->writeLineLoop(_n_lines, _n_plane_sfc, out);
        polygon_tree->writeSubPolygonsAsLineConstraints(_n_lines, _n_plane_sfc-1, out);
        polygon_tree->writeLineConstraints(_n_lines, _n_plane_sfc-1, out);
        polygon_tree->writeStations(pnt_id_offset, _n_plane_sfc-1, out);
        polygon_tree->writeAdditionalPointData(pnt_id_offset, _n_plane_sfc-1, out);
    }

    if (! _keep_preprocessed_geometry) {
        _geo_objs.removeSurfaceVec(_gmsh_geo_name);
        _geo_objs.removePolylineVec(_gmsh_geo_name);
        _geo_objs.removePointVec(_gmsh_geo_name);
        _geo_objs.removeStationVec(gmsh_stations_name);
    }

    return 0;
}

void GMSHInterface::writePoints(std::ostream& out) const
{
    for (auto & gmsh_pnt : _gmsh_pnts) {
        // reverse rotation
        if (gmsh_pnt) {
            double* tmp = _inverse_rot_mat * gmsh_pnt->getCoords();
            (*gmsh_pnt)[0] = tmp[0];
            (*gmsh_pnt)[1] = tmp[1];
            (*gmsh_pnt)[2] = tmp[2];
            delete [] tmp;
            out << *gmsh_pnt << "\n";
        }
    }
}

} // end namespace GMSH
} // end namespace FileIO
