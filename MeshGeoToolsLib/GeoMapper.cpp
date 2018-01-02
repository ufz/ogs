/**
 * \file
 * \author Karsten Rink
 * \date   2012-09-25
 * \brief  Implementation of the GeoMapper class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
}

void GeoMapper::mapOnDEM(std::unique_ptr<GeoLib::Raster const> raster)
{
    std::vector<GeoLib::Point*> const* pnts(_geo_objects.getPointVec(_geo_name));
    if (! pnts) {
        ERR("Geometry \"%s\" does not exist.", _geo_name.c_str());
        return;
    }
    _raster = std::move(raster);

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
        }
        // clean up the copied nodes
        for (std::size_t k(0); k < elem_2d->getNumberOfNodes(); ++k)
        {
            delete elem_2d->getNode(k);
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
            true};
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
            sub_segments.emplace_back(new GeoLib::Point{beg_pnt, 0},
                                      new GeoLib::Point{intersections[0], 0},
                                      true);
    }

    if (intersections.size() == 1 && elem == end_elem)
    {
        // The line segment intersects the element that contains the end
        // point of the line segment. Here the last sub line segment is
        // added.
        if (MathLib::sqrDist(end_pnt, intersections[0]) >
            std::numeric_limits<double>::epsilon())
            sub_segments.emplace_back(new GeoLib::Point{intersections[0], 0},
                                      new GeoLib::Point{end_pnt, 0}, true);
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
        sub_segments.emplace_back(new GeoLib::Point{intersections[0], 0},
                                  new GeoLib::Point{intersections[1], 0}, true);
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
        sub_segments.emplace_back(new GeoLib::Point{beg_pnt, 0}, pnt, true);
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
    GeoLib::LineSegment const& segment)
{
    GeoLib::LineSegment seg_deep_copy(
        new GeoLib::Point(segment.getBeginPoint()),
        new GeoLib::Point(segment.getEndPoint()), true);
    // modify z coordinates such that all surface elements around the line
    // segment are found
    seg_deep_copy.getBeginPoint()[2] = mesh_element_grid.getMinPoint()[2];
    seg_deep_copy.getEndPoint()[2] = mesh_element_grid.getMaxPoint()[2];
    std::array<MathLib::Point3d, 2> const pnts{
        {seg_deep_copy.getBeginPoint(), seg_deep_copy.getEndPoint()}};
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
    for (auto org_line : *org_lines)
    {
        mapPolylineOnSurfaceMesh(*org_line, *org_points, mesh_element_grid);
    }
}

} // end namespace MeshGeoToolsLib
