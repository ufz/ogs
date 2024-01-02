/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Applications/FileIO/Gmsh/GMSHInterface.h"

#include <fstream>
#include <memory>
#include <vector>

#include "Applications/FileIO/Gmsh/GMSHAdaptiveMeshDensity.h"
#include "Applications/FileIO/Gmsh/GMSHFixedMeshDensity.h"
#include "Applications/FileIO/Gmsh/GMSHMeshDensityStrategy.h"
#include "Applications/FileIO/Gmsh/GMSHPoint.h"
#include "Applications/FileIO/Gmsh/GMSHPolygonTree.h"
#include "BaseLib/Algorithm.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polygon.h"
#include "GeoLib/PolygonWithSegmentMarker.h"
#include "GeoLib/PolylineWithSegmentMarker.h"
#include "GeoLib/Station.h"
#include "InfoLib/GitInfo.h"

namespace FileIO
{
namespace GMSH
{
static std::ostream& operator<<(std::ostream& os,
                                std::vector<GMSHPoint*> const& points)
{
    for (auto const* const point : points)
    {
        if (point)
        {
            os << *point << "\n";
        }
    }
    return os;
}

GMSHInterface::GMSHInterface(
    GeoLib::GEOObjects& geo_objs, bool /*include_stations_as_constraints*/,
    GMSH::MeshDensityAlgorithm const mesh_density_algorithm,
    double const pnt_density, double const station_density,
    std::size_t const max_pnts_per_leaf,
    std::vector<std::string> const& selected_geometries, bool const rotate,
    bool const keep_preprocessed_geometry)
    : _n_lines(0),
      _n_plane_sfc(0),
      _geo_objs(geo_objs),
      _selected_geometries(selected_geometries),
      _rotate(rotate),
      _keep_preprocessed_geometry(keep_preprocessed_geometry)
{
    switch (mesh_density_algorithm)
    {
        case GMSH::MeshDensityAlgorithm::FixedMeshDensity:
            _mesh_density_strategy =
                std::make_unique<GMSH::GMSHFixedMeshDensity>(pnt_density);
            break;
        case GMSH::MeshDensityAlgorithm::AdaptiveMeshDensity:
            _mesh_density_strategy =
                std::make_unique<GMSH::GMSHAdaptiveMeshDensity>(
                    pnt_density, station_density, max_pnts_per_leaf);
            break;
    }
}

GMSHInterface::~GMSHInterface()
{
    BaseLib::cleanupVectorElements(_gmsh_pnts);
    for (auto const* polygon_tree : _polygon_tree_list)
    {
        delete polygon_tree;
    }
}

bool GMSHInterface::write()
{
    out << "// GMSH input file created by OpenGeoSys "
        << GitInfoLib::GitInfo::ogs_version;
    out << "\n\n";

    return writeGMSHInputFile(out) <= 0;
}

int GMSHInterface::writeGMSHInputFile(std::ostream& out)
{
    DBUG("GMSHInterface::writeGMSHInputFile(): get data from GEOObjects.");

    if (_selected_geometries.empty())
    {
        return 1;
    }

    // *** get and merge data from _geo_objs
    if (_selected_geometries.size() > 1)
    {
        _gmsh_geo_name = "GMSHGeometry";
        if (_geo_objs.mergeGeometries(_selected_geometries, _gmsh_geo_name))
        {
            return 2;
        }
    }
    else
    {
        _gmsh_geo_name = _selected_geometries[0];
        _keep_preprocessed_geometry = true;
    }

    auto* merged_pnts(const_cast<std::vector<GeoLib::Point*>*>(
        _geo_objs.getPointVec(_gmsh_geo_name)));
    if (!merged_pnts)
    {
        ERR("GMSHInterface::writeGMSHInputFile(): Did not found any points.");
        return 2;
    }

    if (_rotate)
    {
        // Rotate points to the x-y-plane.
        _inverse_rot_mat = GeoLib::rotatePointsToXY(*merged_pnts);
        // Compute inverse rotation matrix to reverse the rotation later on.
        _inverse_rot_mat.transposeInPlace();
    }
    else
    {
        // project data on the x-y-plane
        _inverse_rot_mat = Eigen::Matrix3d::Identity();
        for (auto pnt : *merged_pnts)
        {
            (*pnt)[2] = 0.0;
        }
    }

    std::vector<GeoLib::Polyline*> const* merged_plys(
        _geo_objs.getPolylineVec(_gmsh_geo_name));
    DBUG("GMSHInterface::writeGMSHInputFile(): Obtained data.");

    if (!merged_plys)
    {
        ERR("GMSHInterface::writeGMSHInputFile(): Did not find any polylines.");
        return 2;
    }

    // *** compute and insert all intersection points between polylines
    GeoLib::PointVec& pnt_vec(*const_cast<GeoLib::PointVec*>(
        _geo_objs.getPointVecObj(_gmsh_geo_name)));
    GeoLib::computeAndInsertAllIntersectionPoints(
        pnt_vec, *(const_cast<std::vector<GeoLib::Polyline*>*>(merged_plys)));

    std::vector<GeoLib::Polyline*> polygons;
    // for each closed polyline add a PolygonWithSegmentMarker object into
    // polygons
    for (auto polyline : *merged_plys)
    {
        if (!polyline->isClosed())
        {
            continue;
        }
        polygons.push_back(new GeoLib::PolygonWithSegmentMarker(*polyline));
    }
    if (polygons.empty())
    {
        OGS_FATAL("GMSHInterface::writeGMSHInputFile(): no polygons found.");
    }
    // let the polygon memory management be done by GEOObjects
    _geo_objs.appendPolylineVec(polygons, _gmsh_geo_name);
    // create for each polygon a PolygonTree
    std::transform(
        polygons.begin(), polygons.end(),
        std::back_inserter(_polygon_tree_list),
        [this](auto const& polygon)
        {
            return new GMSH::GMSHPolygonTree(
                dynamic_cast<GeoLib::PolygonWithSegmentMarker*>(polygon),
                nullptr, _geo_objs, _gmsh_geo_name, *_mesh_density_strategy);
        });
    DBUG(
        "GMSHInterface::writeGMSHInputFile(): Computed topological hierarchy - "
        "detected {:d} polygons.",
        _polygon_tree_list.size());
    // compute topological hierarchy of polygons
    GeoLib::createPolygonTrees<GMSH::GMSHPolygonTree>(_polygon_tree_list);
    DBUG(
        "GMSHInterface::writeGMSHInputFile(): Computed topological hierarchy - "
        "calculated {:d} polygon trees.",
        _polygon_tree_list.size());

    // *** Mark in each polygon tree the segments shared by two polygons.
    for (auto* polygon_tree : _polygon_tree_list)
    {
        polygon_tree->markSharedSegments();
    }

    // *** insert stations and polylines (except polygons) in the appropriate
    // object of
    //     class GMSHPolygonTree
    // *** insert stations
    std::vector<GeoLib::Point*> gmsh_stations{};
    for (auto const& geometry_name : _selected_geometries)
    {
        auto const* stations(_geo_objs.getStationVec(geometry_name));
        if (stations)
        {
            for (auto* station : *stations)
            {
                bool found(false);
                for (auto it(_polygon_tree_list.begin());
                     it != _polygon_tree_list.end() && !found;
                     ++it)
                {
                    gmsh_stations.emplace_back(new GeoLib::Station(
                        *static_cast<GeoLib::Station*>(station)));
                    if ((*it)->insertStation(gmsh_stations.back()))
                    {
                        found = true;
                    }
                }
            }
        }
    }
    std::string gmsh_stations_name(_gmsh_geo_name + "-Stations");
    if (!gmsh_stations.empty())
    {
        _geo_objs.addStationVec(std::move(gmsh_stations), gmsh_stations_name);
    }

    // *** insert polylines
    for (auto polyline : *merged_plys)
    {
        if (!polyline->isClosed())
        {
            for (auto* polygon_tree : _polygon_tree_list)
            {
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
    for (std::size_t k(0); k < n_merged_pnts; k++)
    {
        _gmsh_pnts[k] = nullptr;
    }
    for (auto& polygon_tree : _polygon_tree_list)
    {
        polygon_tree->createGMSHPoints(_gmsh_pnts);
    }

    std::stringstream error_messages;
    error_messages.precision(std::numeric_limits<double>::digits10);
    for (std::size_t k = 0; k < _gmsh_pnts.size(); ++k)
    {
        if (_gmsh_pnts[k] == nullptr)
        {
            error_messages
                << "The point at (" << *(*merged_pnts)[k]
                << ") is not part of a polyline, and won't be used in the "
                   "meshing as a constraint. If you want to include it in the "
                   "mesh please create a observation/measurement station for "
                   "the point and include it additional in the meshing "
                   "process.\n";
        }
    }
    auto const error_message = error_messages.str();
    if (!error_message.empty())
    {
        OGS_FATAL("{}", error_message);
    }

    // *** finally write data :-)
    GeoLib::rotatePoints(_inverse_rot_mat, _gmsh_pnts);
    out << _gmsh_pnts;

    std::size_t pnt_id_offset(_gmsh_pnts.size());
    for (auto* polygon_tree : _polygon_tree_list)
    {
        polygon_tree->writeLineLoop(_n_lines, _n_plane_sfc, out,
                                    _write_physical_groups);
        polygon_tree->writeSubPolygonsAsLineConstraints(_n_lines,
                                                        _n_plane_sfc - 1, out);
        polygon_tree->writeLineConstraints(_n_lines, _n_plane_sfc - 1, out);
        polygon_tree->writeStations(pnt_id_offset, _n_plane_sfc - 1, out);
        polygon_tree->writeAdditionalPointData(pnt_id_offset, _n_plane_sfc - 1,
                                               out);
    }

    if (!_keep_preprocessed_geometry)
    {
        _geo_objs.removeSurfaceVec(_gmsh_geo_name);
        _geo_objs.removePolylineVec(_gmsh_geo_name);
        _geo_objs.removePointVec(_gmsh_geo_name);
        _geo_objs.removeStationVec(gmsh_stations_name);
    }

    return 0;
}

}  // end namespace GMSH
}  // end namespace FileIO
