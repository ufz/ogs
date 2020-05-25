/**
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <fstream>
#include <memory>
#include <vector>

#include "BaseLib/Logging.h"

#include "InfoLib/GitInfo.h"
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
GMSHInterface::GMSHInterface(
    GeoLib::GEOObjects& geo_objs, bool /*include_stations_as_constraints*/,
    GMSH::MeshDensityAlgorithm const mesh_density_algorithm,
    double const pnt_density, double const station_density,
    std::size_t const max_pnts_per_leaf,
    std::vector<std::string> const& selected_geometries, bool const rotate,
    bool const keep_preprocessed_geometry)
    : n_lines_(0),
      n_plane_sfc_(0),
      geo_objs_(geo_objs),
      selected_geometries_(selected_geometries),
      rotate_(rotate),
      keep_preprocessed_geometry_(keep_preprocessed_geometry)
{
    switch (mesh_density_algorithm) {
    case GMSH::MeshDensityAlgorithm::FixedMeshDensity:
        mesh_density_strategy_ =
            std::make_unique<GMSH::GMSHFixedMeshDensity>(pnt_density);
        break;
    case GMSH::MeshDensityAlgorithm::AdaptiveMeshDensity:
        mesh_density_strategy_ =
            std::make_unique<GMSH::GMSHAdaptiveMeshDensity>(
                pnt_density, station_density, max_pnts_per_leaf);
        break;
    }
}

GMSHInterface::~GMSHInterface()
{
    for (auto* gmsh_pnt : gmsh_pnts_)
    {
        delete gmsh_pnt;
    }
    for (auto* polygon_tree : polygon_tree_list_)
    {
        delete polygon_tree;
    }
}

bool GMSHInterface::write()
{
    out_ << "// GMSH input file created by OpenGeoSys " << GitInfoLib::GitInfo::ogs_version;
    out_ << "\n\n";

    return writeGMSHInputFile(out_) <= 0;
}

int GMSHInterface::writeGMSHInputFile(std::ostream& out)
{
    DBUG("GMSHInterface::writeGMSHInputFile(): get data from GEOObjects.");

    if (selected_geometries_.empty())
    {
        return 1;
    }

    // *** get and merge data from geo_objs_
    if (selected_geometries_.size() > 1) {
        gmsh_geo_name_ = "GMSHGeometry";
        if (geo_objs_.mergeGeometries(selected_geometries_, gmsh_geo_name_))
        {
            return 2;
        }
    } else {
        gmsh_geo_name_ = selected_geometries_[0];
        keep_preprocessed_geometry_ = true;
    }

    auto* merged_pnts(const_cast<std::vector<GeoLib::Point*>*>(
        geo_objs_.getPointVec(gmsh_geo_name_)));
    if (! merged_pnts) {
        ERR("GMSHInterface::writeGMSHInputFile(): Did not found any points.");
        return 2;
    }

    if (rotate_) {
        // Rotate points to the x-y-plane.
        inverse_rot_mat_ = GeoLib::rotatePointsToXY(*merged_pnts);
        // Compute inverse rotation matrix to reverse the rotation later on.
        inverse_rot_mat_.transposeInPlace();
    } else {
        // project data on the x-y-plane
        inverse_rot_mat_(0,0) = 1.0;
        inverse_rot_mat_(1,1) = 1.0;
        inverse_rot_mat_(2,2) = 1.0;
        for (auto pnt : *merged_pnts)
        {
            (*pnt)[2] = 0.0;
        }
    }

    std::vector<GeoLib::Polyline*> const* merged_plys(
        geo_objs_.getPolylineVec(gmsh_geo_name_));
    DBUG("GMSHInterface::writeGMSHInputFile(): Obtained data.");

    if (!merged_plys) {
        ERR("GMSHInterface::writeGMSHInputFile(): Did not find any polylines.");
        return 2;
    }

    // *** compute and insert all intersection points between polylines
    GeoLib::PointVec& pnt_vec(*const_cast<GeoLib::PointVec*>(
        geo_objs_.getPointVecObj(gmsh_geo_name_)));
    GeoLib::computeAndInsertAllIntersectionPoints(pnt_vec,
        *(const_cast<std::vector<GeoLib::Polyline*>*>(merged_plys)));

    // *** compute topological hierarchy of polygons
    for (auto polyline : *merged_plys) {
        if (!polyline->isClosed()) {
            continue;
        }
        polygon_tree_list_.push_back(new GMSH::GMSHPolygonTree(
            new GeoLib::PolygonWithSegmentMarker(*polyline), nullptr, geo_objs_,
            gmsh_geo_name_, *mesh_density_strategy_));
    }
    DBUG(
        "GMSHInterface::writeGMSHInputFile(): Computed topological hierarchy - "
        "detected {:d} polygons.",
        polygon_tree_list_.size());
    GeoLib::createPolygonTrees<GMSH::GMSHPolygonTree>(polygon_tree_list_);
    DBUG(
        "GMSHInterface::writeGMSHInputFile(): Computed topological hierarchy - "
        "calculated {:d} polygon trees.",
        polygon_tree_list_.size());

    // *** Mark in each polygon tree the segments shared by two polygons.
    for (auto* polygon_tree : polygon_tree_list_)
    {
        polygon_tree->markSharedSegments();
    }

    // *** insert stations and polylines (except polygons) in the appropriate object of
    //     class GMSHPolygonTree
    // *** insert stations
    auto gmsh_stations = std::make_unique<std::vector<GeoLib::Point*>>();
    for (auto const& geometry_name : selected_geometries_) {
        auto const* stations(geo_objs_.getStationVec(geometry_name));
        if (stations) {
            for (auto * station : *stations) {
                bool found(false);
                for (auto it(polygon_tree_list_.begin());
                    it != polygon_tree_list_.end() && !found; ++it) {
                    gmsh_stations->emplace_back(new GeoLib::Station(
                        *static_cast<GeoLib::Station*>(station)));
                    if ((*it)->insertStation(gmsh_stations->back())) {
                        found = true;
                    }
                }
            }
        }
    }
    std::string gmsh_stations_name(gmsh_geo_name_+"-Stations");
    if (! gmsh_stations->empty()) {
        geo_objs_.addStationVec(std::move(gmsh_stations), gmsh_stations_name);
    }

    // *** insert polylines
    for (auto polyline : *merged_plys) {
        if (!polyline->isClosed()) {
            for (auto * polygon_tree : polygon_tree_list_) {
                auto polyline_with_segment_marker =
                    new GeoLib::PolylineWithSegmentMarker(*polyline);
                polygon_tree->insertPolyline(polyline_with_segment_marker);
            }
        }
    }

    // *** init mesh density strategies
    for (auto& polygon_tree : polygon_tree_list_)
    {
        polygon_tree->initMeshDensityStrategy();
    }

    // *** create GMSH data structures
    const std::size_t n_merged_pnts(merged_pnts->size());
    gmsh_pnts_.resize(n_merged_pnts);
    for (std::size_t k(0); k<n_merged_pnts; k++) {
        gmsh_pnts_[k] = nullptr;
    }
    for (auto& polygon_tree : polygon_tree_list_)
    {
        polygon_tree->createGMSHPoints(gmsh_pnts_);
    }

    // *** finally write data :-)
    writePoints(out);
    std::size_t pnt_id_offset(gmsh_pnts_.size());
    for (auto* polygon_tree : polygon_tree_list_)
    {
        polygon_tree->writeLineLoop(n_lines_, n_plane_sfc_, out);
        polygon_tree->writeSubPolygonsAsLineConstraints(n_lines_, n_plane_sfc_-1, out);
        polygon_tree->writeLineConstraints(n_lines_, n_plane_sfc_-1, out);
        polygon_tree->writeStations(pnt_id_offset, n_plane_sfc_-1, out);
        polygon_tree->writeAdditionalPointData(pnt_id_offset, n_plane_sfc_-1, out);
    }

    if (! keep_preprocessed_geometry_) {
        geo_objs_.removeSurfaceVec(gmsh_geo_name_);
        geo_objs_.removePolylineVec(gmsh_geo_name_);
        geo_objs_.removePointVec(gmsh_geo_name_);
        geo_objs_.removeStationVec(gmsh_stations_name);
    }

    return 0;
}

void GMSHInterface::writePoints(std::ostream& out) const
{
    for (auto & gmsh_pnt : gmsh_pnts_) {
        // reverse rotation
        if (gmsh_pnt) {
            double* tmp = inverse_rot_mat_ * gmsh_pnt->getCoords();
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
