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
#include <sstream>
#include <numeric>

#include <logog/include/logog.hpp>

#include "BaseLib/makeVectorUnique.h"
#include "BaseLib/Error.h"

#include "GeoLib/AABB.h"
#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Raster.h"
#include "GeoLib/StationBorehole.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/FaceRule.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshSurfaceExtraction.h"
#include "MeshLib/MeshEditing/projectMeshOntoPlane.h"
#include "MeshLib/MeshSearch/MeshElementGrid.h"

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

/// Find the 2d-element within the \c elements that contains the given point \c p.
///
/// The algorithm projects every element of the elements vector and the point
/// \c p orthogonal to the \f$x\f$-\f$y\f$ plane. In the \f$x\f$-\f$y\f$ plane
/// it is checked if the projected point is in the projected element.
static MeshLib::Element const* findElementContainingPointXY(
    std::vector<MeshLib::Element const*> const& elements,
    MathLib::Point3d const& p)
{
    for (auto const elem : elements) {
        std::unique_ptr<MeshLib::Element> elem_2d(elem->clone());
        // reset/copy the nodes
        for (std::size_t k(0); k<elem_2d->getNumberOfNodes(); ++k) {
            elem_2d->setNode(k, new MeshLib::Node(*elem_2d->getNode(k)));
        }
        // project to xy
        for (std::size_t k(0); k<elem_2d->getNumberOfNodes(); ++k) {
            (*const_cast<MeshLib::Node*>(elem_2d->getNode(k)))[2] = 0.0;
        }
        if (elem_2d->isPntInElement(MathLib::Point3d{ {{p[0], p[1], 0.0}} })) {
            // clean up the copied nodes
            for (std::size_t k(0); k<elem_2d->getNumberOfNodes(); ++k) {
                delete elem_2d->getNode(k);
            }
            return elem;
        } else {
            // clean up the copied nodes
            for (std::size_t k(0); k<elem_2d->getNumberOfNodes(); ++k) {
                delete elem_2d->getNode(k);
            }
        }
    }
    return nullptr;
}

static std::vector<MathLib::Point3d> computeElementSegmentIntersections(
    MeshLib::Element const& elem, GeoLib::LineSegment const& segment)
{
    std::vector<MathLib::Point3d> element_intersections;
    for (std::size_t k(0); k < elem.getNumberOfEdges(); ++k)
    {
        auto const edge =
            std::unique_ptr<MeshLib::Element const>(elem.getEdge(k));
        GeoLib::LineSegment elem_segment{
            new GeoLib::Point(*dynamic_cast<MathLib::Point3d*>(
                const_cast<MeshLib::Node*>(edge->getNode(0))), 0),
            new GeoLib::Point(*dynamic_cast<MathLib::Point3d*>(
                const_cast<MeshLib::Node*>(edge->getNode(1))), 0),
            false};
        std::vector<MathLib::Point3d> const intersections(
            GeoLib::lineSegmentIntersect2d(segment, elem_segment));
        for (auto const& p : intersections)
            element_intersections.push_back(std::move(p));
    }
    return element_intersections;
}

static std::vector<GeoLib::LineSegment> createSubSegmentsForElement(
    std::vector<MathLib::Point3d> const& intersections,
    MeshLib::Element const* const beg_elem,
    MeshLib::Element const* const end_elem, MathLib::Point3d const& beg_pnt,
    MathLib::Point3d const& end_pnt, MeshLib::Element const* const elem)
{
    std::vector<GeoLib::LineSegment> sub_segments;
    if (intersections.size() > 2) {
        std::stringstream out;
        out << "element with id " << elem->getID() << " and seg "
                  << " intersecting at more than two edges\n";
        for (std::size_t k(0); k<intersections.size(); ++k)
            out << k << " " << intersections[k]
                      << "\n";
        out << "Could not map segment on element. Aborting.\n";
        OGS_FATAL("%s", out.str().c_str());
    }

    if (intersections.size() == 1 && elem == beg_elem)
    {
        // The line segment intersects the element that contains the begin
        // point of the line segment. Here the first sub line segment is
        // added.
        if (MathLib::sqrDist(beg_pnt, intersections[0]) >
            std::numeric_limits<double>::epsilon())
            sub_segments.emplace_back(GeoLib::LineSegment{
                new GeoLib::Point{beg_pnt, 0},
                new GeoLib::Point{intersections[0], 0}, true});
    }

    if (intersections.size() == 1 && elem == end_elem)
    {
        // The line segment intersects the element that contains the end
        // point of the line segment. Here the last sub line segment is
        // added.
        if (MathLib::sqrDist(end_pnt, intersections[0]) >
            std::numeric_limits<double>::epsilon())
            sub_segments.emplace_back(
                GeoLib::LineSegment{new GeoLib::Point{intersections[0], 0},
                                    new GeoLib::Point{end_pnt, 0}, true});
    }

    if (intersections.size() == 1 && (elem != beg_elem && elem != end_elem))
    {
        // Since the line segment enters and leaves the element in the same
        // point there isn't any need to insert a new sub line segment.
        return sub_segments;
    }

    // create sub segment for the current element
    if (intersections.size() == 2)
    {
        sub_segments.emplace_back(
            GeoLib::LineSegment{new GeoLib::Point{intersections[0], 0},
                                new GeoLib::Point{intersections[1], 0}, true});
    }
    return sub_segments;
}

static std::vector<GeoLib::LineSegment> mapLineSegment(
    GeoLib::LineSegment const& segment,
    std::vector<MeshLib::Element const*> const& surface_elements,
    MeshLib::Element const* const beg_elem,
    MeshLib::Element const* const end_elem)
{
    std::vector<GeoLib::LineSegment> sub_segments;
    MathLib::Point3d const& beg_pnt(segment.getBeginPoint());
    MathLib::Point3d const& end_pnt(segment.getEndPoint());

    for (auto const elem : surface_elements) {
        // compute element-segment-intersections (2d in x-y-plane)
        std::vector<MathLib::Point3d> element_intersections(
            computeElementSegmentIntersections(*elem, segment));
        if (element_intersections.empty())
            continue;

        BaseLib::makeVectorUnique(element_intersections);

        std::vector<GeoLib::LineSegment> sub_seg_elem(
            createSubSegmentsForElement(element_intersections, beg_elem,
                                        end_elem, beg_pnt, end_pnt, elem));
        sub_segments.insert(sub_segments.end(), sub_seg_elem.begin(),
                            sub_seg_elem.end());
    }

    // beg_elem == nullptr means there isn't any element corresponding to the
    // beg_pnt and as a consequence the above algorithm doesn't insert a sub
    // segment
    if (beg_elem == nullptr)
    {
        auto min_dist_segment = std::min_element(
            sub_segments.begin(), sub_segments.end(),
            [&beg_pnt](GeoLib::LineSegment const& seg0,
                       GeoLib::LineSegment const& seg1)
            {
                // min dist for segment 0
                const double d0(
                    std::min(MathLib::sqrDist(beg_pnt, seg0.getBeginPoint()),
                             MathLib::sqrDist(beg_pnt, seg0.getEndPoint())));
                // min dist for segment 1
                const double d1(
                    std::min(MathLib::sqrDist(beg_pnt, seg1.getBeginPoint()),
                             MathLib::sqrDist(beg_pnt, seg1.getEndPoint())));
                return d0 < d1;
            });
        GeoLib::Point * pnt{
            MathLib::sqrDist(beg_pnt, min_dist_segment->getBeginPoint()) <
                    MathLib::sqrDist(beg_pnt, min_dist_segment->getEndPoint())
                ? new GeoLib::Point{min_dist_segment->getBeginPoint()}
                : new GeoLib::Point{min_dist_segment->getEndPoint()}};
        sub_segments.emplace_back(
            GeoLib::LineSegment{new GeoLib::Point{beg_pnt, 0}, pnt, true});
    }
    // sort all sub segments for the given segment (beg_pnt, end_pnt)
    GeoLib::sortSegments(beg_pnt, sub_segments);

    sub_segments.erase(std::unique(sub_segments.begin(), sub_segments.end()),
                       sub_segments.end());

    return sub_segments;
}

static void mapPointOnSurfaceElement(MeshLib::Element const& elem,
                                     MathLib::Point3d& q)
{
    // create plane equation: n*p = d
    MathLib::Vector3 const& p(*(elem.getNode(0)));
    MathLib::Vector3 const n(MeshLib::FaceRule::getSurfaceNormal(&elem));
    if (n[2] == 0.0) { // vertical plane, z coordinate is arbitrary
        q[2] = p[2];
    } else {
        double const d(MathLib::scalarProduct(n, p));
        q[2] = (d - n[0]*q[0] - n[1]*q[1])/n[2];
    }
}

static std::vector<MeshLib::Element const*>
getCandidateElementsForLineSegmentIntersection(
    MeshLib::MeshElementGrid const& mesh_element_grid,
    GeoLib::LineSegment segment)
{
    // modify z coordinates such that all surface elements around the line
    // segment are found
    segment.getBeginPoint()[2] = mesh_element_grid.getMinPoint()[2];
    segment.getEndPoint()[2] = mesh_element_grid.getMaxPoint()[2];
    std::array<MathLib::Point3d, 2> const pnts{
        {segment.getBeginPoint(), segment.getEndPoint()}};
    GeoLib::AABB aabb(pnts.cbegin(), pnts.cend());

    auto candidate_elements = mesh_element_grid.getElementsInVolume(
        aabb.getMinPoint(), aabb.getMaxPoint());

    // make candidate elements unique
    BaseLib::makeVectorUnique(candidate_elements);

    return candidate_elements;
}

static bool snapPointToElementNode(MathLib::Point3d& p,
                                   MeshLib::Element const& elem, double rel_eps)
{
    // values will be initialized within computeSqrNodeDistanceRange
    double sqr_min, sqr_max;
    elem.computeSqrNodeDistanceRange(sqr_min, sqr_max);

    double const sqr_eps(rel_eps*rel_eps * sqr_min);
    for (std::size_t k(0); k<elem.getNumberOfNodes(); ++k) {
        auto const& node(*elem.getNode(k));
        double const sqr_dist_2d(MathLib::sqrDist2d(p, node));
        if (sqr_dist_2d < sqr_eps) {
#ifdef DEBUG_GEOMAPPER
            std::stringstream out;
            out.precision(std::numeric_limits<double>::digits10);
            out << "Segment point snapped from " << p;
#endif
            p = node;
#ifdef DEBUG_GEOMAPPER
            out << "to " << p;
            DBUG("%s", out.str().c_str());
#endif
            return true;
        }
    }
    return false;
}

static void insertSubSegments(
    GeoLib::Polyline& ply, GeoLib::PointVec& points,
    GeoLib::Polyline::SegmentIterator& segment_it,
    std::vector<GeoLib::LineSegment> const& sub_segments)
{
    std::size_t const j(segment_it.getSegmentNumber());
    std::size_t new_pnts_cnt(0);
    for (auto const& segment : sub_segments)
    {
        auto const begin_id(points.push_back(
            new GeoLib::Point(segment.getBeginPoint(), points.size())));
        if (ply.insertPoint(j + new_pnts_cnt + 1, begin_id))
            new_pnts_cnt++;
        auto const end_id(points.push_back(
            new GeoLib::Point(segment.getEndPoint(), points.size())));
        if (ply.insertPoint(j + new_pnts_cnt + 1, end_id))
            new_pnts_cnt++;
    }
    segment_it += new_pnts_cnt;
}

static void mapPolylineOnSurfaceMesh(
    GeoLib::Polyline& ply,
    GeoLib::PointVec& orig_points,
    MeshLib::MeshElementGrid const& mesh_element_grid)
{
    // for each segment ...
    for (auto segment_it(ply.begin()); segment_it != ply.end(); ++segment_it)
    {
        auto candidate_elements(getCandidateElementsForLineSegmentIntersection(
            mesh_element_grid, *segment_it));

        auto mapPoint = [&candidate_elements](MathLib::Point3d& p) {
            auto const* elem(
                findElementContainingPointXY(candidate_elements, p));
            if (elem)
            {
                if (!snapPointToElementNode(p, *elem, 1e-3))
                    mapPointOnSurfaceElement(*elem, p);
            }
            return elem;
        };

        // map segment begin and end point
        auto const* beg_elem(mapPoint((*segment_it).getBeginPoint()));
        auto const* end_elem(mapPoint((*segment_it).getEndPoint()));

        // Since the mapping of the segment begin and end points the coordinates
        // changed. The internal data structures of PointVec are possibly
        // invalid and hence it is necessary to re-create them.
        orig_points.resetInternalDataStructures();

        if (beg_elem == end_elem) {
            // TODO: handle cases: beg_elem == end_elem == nullptr
            // There are further checks necessary to determine which case we are
            // in:
            // 1. beg_elem == end_elem and the segment intersects elements
            // 2. beg_elem == end_elem and the segment does not intersect any
            // element, i.e., the segment is located outside of the mesh area
            //
            // Case 1 needs additional work.
            continue;
        }

        // map the line segment (and if necessary for the mapping partition it)
        std::vector<GeoLib::LineSegment> sub_segments(mapLineSegment(
            *segment_it, candidate_elements, beg_elem, end_elem));

        if (sub_segments.empty())
            continue;

        // The case sub_segment.size() == 1 is already handled above.

        if (sub_segments.size() > 1)
        {
            insertSubSegments(ply, orig_points, segment_it, sub_segments);
        }
    }
}

void GeoMapper::advancedMapOnMesh(MeshLib::Mesh const& mesh)
{
    // 1. extract surface
    if (_surface_mesh)
        delete _surface_mesh;
    if (mesh.getDimension()<3) {
        _surface_mesh = new MeshLib::Mesh(mesh);
    } else {
        const MathLib::Vector3 dir(0,0,-1);
        _surface_mesh =
            MeshLib::MeshSurfaceExtraction::getMeshSurface(mesh, dir, 90+1e-6);
    }

    // 2. compute mesh grid for surface
    MeshLib::MeshElementGrid const mesh_element_grid(*_surface_mesh);

    // 3. map each polyline
    auto org_lines(_geo_objects.getPolylineVec(_geo_name));
    auto org_points(_geo_objects.getPointVecObj(_geo_name));
    for (std::size_t k(0); k<org_lines->size(); ++k) {
        mapPolylineOnSurfaceMesh(*((*org_lines)[k]), *org_points,
                                 mesh_element_grid);
    }
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
