/**
 * \file
 * \author Karsten Rink
 * \date   2012-09-25
 * \brief  Implementation of the GeoMapper class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GeoMapper.h"

#include <algorithm>
#include <numeric>

#include <logog/include/logog.hpp>

#include "GeoLib/AABB.h"
#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Raster.h"
#include "GeoLib/StationBorehole.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshSurfaceExtraction.h"
#include "MeshLib/MeshEditing/projectMeshOntoPlane.h"
#include "MeshLib/IO/readMeshFromFile.h"

namespace MeshGeoToolsLib {

GeoMapper::GeoMapper(GeoLib::GEOObjects &geo_objects, const std::string &geo_name)
    : _geo_objects(geo_objects), _geo_name(const_cast<std::string&>(geo_name)),
    _surface_mesh(nullptr), _grid(nullptr), _raster(nullptr)
{
}

GeoMapper::~GeoMapper()
{
    delete _surface_mesh;
    delete _raster;
}

void GeoMapper::mapOnDEM(GeoLib::Raster *const raster)
{
    std::vector<GeoLib::Point*> const* pnts(_geo_objects.getPointVec(_geo_name));
    if (! pnts) {
        ERR("Geometry \"%s\" does not exist.", _geo_name.c_str());
        return;
    }
    _raster = raster;

    if (GeoLib::isStation((*pnts)[0])) {
        mapStationData(*pnts);
    } else {
        mapPointDataToDEM(*pnts);
    }
}

void GeoMapper::mapOnMesh(MeshLib::Mesh const*const mesh)
{
    std::vector<GeoLib::Point*> const* pnts(_geo_objects.getPointVec(_geo_name));
    if (! pnts) {
        ERR("Geometry \"%s\" does not exist.", _geo_name.c_str());
        return;
    }

    // the variable _surface_mesh is reused below, so first the existing
    // _surface_mesh has to be cleaned up
    if (_surface_mesh)
        delete _surface_mesh;

    if (mesh->getDimension()<3)
        _surface_mesh = new MeshLib::Mesh(*mesh);
    else
    {
        const MathLib::Vector3 dir(0,0,-1);
        _surface_mesh = MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh, dir, 90);
    }

    // init grid
    MathLib::Point3d origin(std::array<double,3>{{0,0,0}});
    MathLib::Vector3 normal(0,0,-1);
    std::vector<MeshLib::Node> flat_nodes;
    flat_nodes.reserve(_surface_mesh->getNumberOfNodes());
    // copy nodes and project the copied nodes to the x-y-plane, i.e. set
    // z-coordinate to zero
    for (auto n_ptr : _surface_mesh->getNodes()) {
        flat_nodes.emplace_back(*n_ptr);
        flat_nodes.back()[2] = 0.0;
    }
    _grid = new GeoLib::Grid<MeshLib::Node>(flat_nodes.cbegin(), flat_nodes.cend());

    if (GeoLib::isStation((*pnts)[0])) {
        mapStationData(*pnts);
    } else {
        mapPointDataToMeshSurface(*pnts);
    }

    delete _grid;
}

void GeoMapper::mapToConstantValue(double value)
{
    std::vector<GeoLib::Point*> const* points (this->_geo_objects.getPointVec(this->_geo_name));
    if (points == nullptr)
    {
        ERR ("Geometry \"%s\" not found.", this->_geo_name.c_str());
        return;
    }
    std::for_each(points->begin(), points->end(), [value](GeoLib::Point* pnt){ (*pnt)[2] = value; });
}

void GeoMapper::mapStationData(std::vector<GeoLib::Point*> const& points)
{
    double min_val(0), max_val(0);
    if (_surface_mesh)
    {
        GeoLib::AABB bounding_box(
            _surface_mesh->getNodes().begin(), _surface_mesh->getNodes().end());
        min_val = bounding_box.getMinPoint()[2];
        max_val = bounding_box.getMaxPoint()[2];
    }

    for (auto * pnt : points)
    {
        double offset =
            (_grid)
                ? (getMeshElevation((*pnt)[0], (*pnt)[1], min_val, max_val) -
                   (*pnt)[2])
                : getDemElevation(*pnt);

        if (!GeoLib::isBorehole(pnt))
            continue;
        auto const& layers = static_cast<GeoLib::StationBorehole*>(pnt)->getProfile();
        for (auto * layer_pnt : layers) {
            (*layer_pnt)[2] = (*layer_pnt)[2] + offset;
        }
    }
}

void GeoMapper::mapPointDataToDEM(std::vector<GeoLib::Point*> const& points)
{
    for (auto * pnt : points)
    {
        GeoLib::Point &p(*pnt);
        p[2] = getDemElevation(p);
    }
}

void GeoMapper::mapPointDataToMeshSurface(std::vector<GeoLib::Point*> const& pnts)
{
    GeoLib::AABB const aabb(
        _surface_mesh->getNodes().cbegin(), _surface_mesh->getNodes().cend());
    double const min_val(aabb.getMinPoint()[2]);
    double const max_val(aabb.getMaxPoint()[2]);

    for (auto * pnt : pnts) {
        // check if pnt is inside of the bounding box of the _surface_mesh
        // projected onto the y-x plane
        GeoLib::Point &p(*pnt);
        if (p[0] < aabb.getMinPoint()[0] || aabb.getMaxPoint()[0] < p[0])
            continue;
        if (p[1] < aabb.getMinPoint()[1] || aabb.getMaxPoint()[1] < p[1])
            continue;

        p[2] = getMeshElevation(p[0], p[1], min_val, max_val);
    }
}

float GeoMapper::getDemElevation(GeoLib::Point const& pnt) const
{
    double const elevation (_raster->getValueAtPoint(pnt));
    if (std::abs(elevation-_raster->getHeader().no_data) < std::numeric_limits<double>::epsilon())
        return 0.0;
    return static_cast<float>(elevation);
}

double GeoMapper::getMeshElevation(
    double x, double y, double min_val, double max_val) const
{
    const MeshLib::Node* pnt =
        _grid->getNearestPoint(MathLib::Point3d{{{x, y, 0}}});
    const std::vector<MeshLib::Element*> elements(
        _surface_mesh->getNode(pnt->getID())->getElements());
    std::unique_ptr<GeoLib::Point> intersection;

    for (auto const & element : elements)
    {
        if (intersection == nullptr &&
            element->getGeomType() != MeshLib::MeshElemType::LINE)
            intersection = GeoLib::triangleLineIntersection(
                *element->getNode(0), *element->getNode(1),
                *element->getNode(2), GeoLib::Point(x, y, max_val),
                GeoLib::Point(x, y, min_val));

        if (intersection == nullptr &&
            element->getGeomType() == MeshLib::MeshElemType::QUAD)
            intersection = GeoLib::triangleLineIntersection(
                *element->getNode(0), *element->getNode(2),
                *element->getNode(3), GeoLib::Point(x, y, max_val),
                GeoLib::Point(x, y, min_val));
    }
    if (intersection)
        return (*intersection)[2];
    // if something goes wrong, simply take the elevation of the nearest mesh node
    return (*(_surface_mesh->getNode(pnt->getID())))[2];
}

std::unique_ptr<std::vector<GeoLib::Polyline*>> copyPolylinesVector(
    std::vector<GeoLib::Polyline*> const& polylines,
    std::vector<GeoLib::Point*> const& points)
{
    std::size_t nLines = polylines.size();
    auto new_lines = std::unique_ptr<std::vector<GeoLib::Polyline*>>(
        new std::vector<GeoLib::Polyline*>(nLines));

    for (std::size_t i=0; i<nLines; ++i)
    {
        (*new_lines)[i] = new GeoLib::Polyline(points);
        std::size_t nLinePnts (polylines[i]->getNumberOfPoints());
        for (std::size_t j=0; j<nLinePnts; ++j)
            (*new_lines)[i]->addPoint(polylines[i]->getPointID(j));
    }
    return new_lines;
}


void GeoMapper::advancedMapOnMesh(
    MeshLib::Mesh const* mesh, std::string const& new_geo_name)
{
    const std::vector<GeoLib::Point*> *points(_geo_objects.getPointVec(_geo_name));
    const std::vector<GeoLib::Polyline*> *org_lines(_geo_objects.getPolylineVec(_geo_name));

    const GeoLib::AABB aabb(points->begin(), points->end());
    const double eps = sqrt(std::numeric_limits<float>::epsilon()) *
                       sqrt(MathLib::sqrDist(aabb.getMinPoint(),aabb.getMaxPoint())) ;

    // copy geometry (and set z=0 for all points)
    auto new_points = std::unique_ptr<std::vector<GeoLib::Point*>>(
        new std::vector<GeoLib::Point*>);
    new_points->reserve(points->size());
    std::transform(points->cbegin(), points->cend(), std::back_inserter(*new_points),
        [](GeoLib::Point* p) { return new GeoLib::Point((*p)[0],(*p)[1],0.0); });

    auto new_lines = copyPolylinesVector(*_geo_objects.getPolylineVec(_geo_name),
                                         *new_points);

    GeoLib::Grid<GeoLib::Point> grid(new_points->begin(), new_points->end());
    double max_segment_length(getMaxSegmentLength(*new_lines));
    // squared so it can be compared to the squared distances calculated later
    max_segment_length *= max_segment_length;

    const unsigned nMeshNodes ( mesh->getNumberOfNodes() );
    // index of closest geo point for each mesh node in (x,y)-plane
    std::vector<int> closest_geo_point(nMeshNodes, -1);
    // distance between geo points and mesh nodes in (x,y)-plane
    std::vector<double> dist(nMeshNodes);
    auto zero_coords = GeoLib::Point{};  // All coordinates zero.
    for (std::size_t i=0; i<nMeshNodes; ++i) {
        zero_coords[0] = (*mesh->getNode(i))[0];
        zero_coords[1] = (*mesh->getNode(i))[1];
        GeoLib::Point* pnt = grid.getNearestPoint(zero_coords);
        dist[i] = MathLib::sqrDist(*pnt, zero_coords);
        closest_geo_point[i] = (dist[i]<=max_segment_length) ? pnt->getID() : -1;
    }

    // store for each point the line segment to which it was added.
    const std::size_t nLines(new_lines->size());
    std::vector< std::vector<unsigned> > line_segment_map(nLines);
    for (std::size_t i=0; i<nLines; ++i)
    {
        line_segment_map[i] = std::vector<unsigned>((*new_lines)[i]->getNumberOfPoints(),0);
        std::iota(line_segment_map[i].begin(), line_segment_map[i].end(), 0);
    }

    for (std::size_t i=0; i<nMeshNodes; ++i)
    {
        // if mesh node too far away or exactly at point position
        if (closest_geo_point[i] == -1 || dist[i] < eps) continue;

        const MeshLib::Node* node (mesh->getNode(i));
        for (std::size_t l=0; l<nLines; ++l)
        {
            // find relevant polylines
            if (!(*org_lines)[l]->isPointIDInPolyline(closest_geo_point[i])) continue;

            // find point position of closest geo point in original polyline
            GeoLib::Polyline* ply ((*org_lines)[l]);
            std::size_t nLinePnts ( ply->getNumberOfPoints() );
            std::size_t node_index_in_ply (0);
            for (node_index_in_ply=0; node_index_in_ply<nLinePnts; ++node_index_in_ply)
                if (ply->getPoint(node_index_in_ply) == (*points)[closest_geo_point[i]])
                    break;
            const GeoLib::Point* geo_point (ply->getPoint(node_index_in_ply));

            // check if line segments connected to closest geo point intersect connected elements of current node
            const std::vector<MeshLib::Element*>& elements (node->getElements());
            const std::size_t nElems = elements.size();
            for (std::size_t e=0; e<nElems; ++e)
            {
                const unsigned nEdges (elements[e]->getNumberOfEdges());
                unsigned intersection_count (0);

                for (unsigned n=0; n<nEdges; ++n)
                {
                    if (intersection_count>1) break; //already two intersections

                    const MeshLib::Element* line = elements[e]->getEdge(n);
                    unsigned index_offset(0); // default: add to first line segment
                    GeoLib::Point* intersection (nullptr);
                    if (node_index_in_ply>0) // test line segment before closest point
                        intersection = calcIntersection(line->getNode(0), line->getNode(1), geo_point, ply->getPoint(node_index_in_ply-1));
                    if (intersection == nullptr && node_index_in_ply<(nLinePnts-1)) // test line segment after closest point
                    {
                        intersection = calcIntersection(line->getNode(0), line->getNode(1), geo_point, ply->getPoint(node_index_in_ply+1));
                        index_offset = 1; // add to second segment
                    }
                    if (intersection)
                    {
                        intersection_count++;
                        unsigned start_point_idx = static_cast<unsigned>(std::distance(line_segment_map[l].begin(), std::find_if(line_segment_map[l].begin(), line_segment_map[l].end(), [&node_index_in_ply, &index_offset](unsigned a){return a==node_index_in_ply+index_offset-1;})));
                        unsigned end_point_idx   = static_cast<unsigned>(std::distance(line_segment_map[l].begin(), std::find_if(line_segment_map[l].begin(), line_segment_map[l].end(), [&node_index_in_ply, &index_offset](unsigned a){return a==node_index_in_ply+index_offset;})));
                        std::size_t pos = getPointPosInLine((*new_lines)[l], start_point_idx, end_point_idx, intersection, eps);

                        if (pos)
                        {
                            const std::size_t pnt_pos (new_points->size());
                            new_points->push_back(intersection);
                            (*new_lines)[l]->insertPoint(pos, pnt_pos);
                            line_segment_map[l].insert(line_segment_map[l].begin()+pos, node_index_in_ply+index_offset-1);
                        }
                    }
                }
            }
        }
    }

    _geo_objects.addPointVec(std::move(new_points), const_cast<std::string&>(new_geo_name));
    std::vector<std::size_t> pnt_id_map = this->_geo_objects.getPointVecObj(new_geo_name)->getIDMap();
    for (auto & new_line : *new_lines)
        new_line->updatePointIDs(pnt_id_map);
    _geo_objects.addPolylineVec(std::move(new_lines), new_geo_name);

    // map new geometry incl. additional point using the normal mapping method
    this->_geo_name = new_geo_name;
    this->mapOnMesh(mesh);
}

GeoLib::Point* GeoMapper::calcIntersection(MathLib::Point3d const*const p1, MathLib::Point3d const*const p2, GeoLib::Point const*const q1, GeoLib::Point const*const q2) const
{
    const double x1 = (*p1)[0], x2 = (*p2)[0], x3 = (*q1)[0], x4 = (*q2)[0];
    const double y1 = (*p1)[1], y2 = (*p2)[1], y3 = (*q1)[1], y4 = (*q2)[1];

    const double det = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if (std::abs(det) < std::numeric_limits<double>::epsilon()) return nullptr;

    const double pre  = (x1*y2 - y1*x2);
    const double post = (x3*y4 - y3*x4);
    const double x = ( pre * (x3 - x4) - (x1 - x2) * post ) / det;
    const double y = ( pre * (y3 - y4) - (y1 - y2) * post ) / det;

    // Check if the x and y coordinates are within both line segments
    if (isPntInBoundingBox(x1,y1,x2,y2,x,y) && isPntInBoundingBox(x3,y3,x4,y4,x,y))
        return new GeoLib::Point(x, y, 0);
    return nullptr;
}

unsigned GeoMapper::getPointPosInLine(GeoLib::Polyline const*const line, unsigned start, unsigned end, GeoLib::Point const*const point, double eps) const
{
    const double* first = line->getPoint(start)->getCoords();
    const double* pnt   = point->getCoords();
    const double max_dist = MathLib::sqrDist(first, pnt);

    // if point is at start or end of line segment
    if (max_dist<eps && MathLib::sqrDist(pnt, line->getPoint(end)->getCoords())) return 0;

    for (std::size_t i=start+1; i<end; ++i)
    {
        const double* current = (*line->getPoint(i)).getCoords();
        if (MathLib::sqrDist(pnt, current) < eps) return 0;
        if (MathLib::sqrDist(first, current) > max_dist) return i;
    }
    return end; // last point of segment
}

bool GeoMapper::isPntInBoundingBox(double ax, double ay, double bx, double by, double px, double py) const
{
    if ( px < (std::min(ax, bx)-std::numeric_limits<double>::epsilon()) || px > (std::max(ax, bx)+std::numeric_limits<double>::epsilon()) ||
         py < (std::min(ay, by)-std::numeric_limits<double>::epsilon()) || py > (std::max(ay, by)+std::numeric_limits<double>::epsilon()) )
         return false;
    return true;
}

double GeoMapper::getMaxSegmentLength(const std::vector<GeoLib::Polyline*> &lines) const
{
    double max_segment_length (0);
    const std::size_t nPlys ( lines.size() );
    for (std::size_t i=0; i<nPlys; ++i)
    {
        const GeoLib::Polyline* line = lines[i];
        const std::size_t nPlyPoints = line->getNumberOfPoints();
        for (std::size_t j=1; j<nPlyPoints; ++j)
        {
            const double dist (line->getLength(j)-line->getLength(j-1));
            if (dist>max_segment_length)
                max_segment_length=dist;
        }
    }
    return max_segment_length;
}

} // end namespace MeshGeoToolsLib

