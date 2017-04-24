/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GMSHPolygonTree.h"

#include "GMSHFixedMeshDensity.h"
#include "GMSHAdaptiveMeshDensity.h"

#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polygon.h"
#include "GeoLib/PolylineWithSegmentMarker.h"
#include "GeoLib/PolygonWithSegmentMarker.h"

namespace FileIO
{
namespace GMSH
{

GMSHPolygonTree::GMSHPolygonTree(GeoLib::PolygonWithSegmentMarker* polygon,
                GMSHPolygonTree* parent,
                GeoLib::GEOObjects &geo_objs, std::string const& geo_name,
                GMSHMeshDensityStrategy * mesh_density_strategy) :
    GeoLib::SimplePolygonTree(polygon, parent), _geo_objs(geo_objs), _geo_name(geo_name),
    _mesh_density_strategy(mesh_density_strategy)
{}

GMSHPolygonTree::~GMSHPolygonTree()
{
    // the polylines are processed also by the children, but the root is
    // responsible to cleanup up
    if (_parent == nullptr) { // root
        for (auto * polyline : _plys)
            delete polyline;
    }
    // member of GeoLib::SimplePolygonTree, but the ownership is not transmitted
    delete _node_polygon;
}

void GMSHPolygonTree::markSharedSegments()
{
    if (_children.empty())
        return;

    if (_parent == nullptr)
        return;

    for (auto& it : _children)
    {
        std::size_t const n_pnts(it->getPolygon()->getNumberOfPoints());
        for (std::size_t k(1); k<n_pnts; k++) {
            if (GeoLib::containsEdge(*(_parent->getPolygon()),
                _node_polygon->getPointID(k-1),
                _node_polygon->getPointID(k)))
            static_cast<GeoLib::PolygonWithSegmentMarker*>(_node_polygon)->markSegment(k, true);
        }
    }
}

bool GMSHPolygonTree::insertStation(GeoLib::Point const* station)
{
    if (_node_polygon->isPntInPolygon(*station)) {
        // try to insert station into the child nodes
        for (std::list<SimplePolygonTree*>::const_iterator it (_children.begin());
             it != _children.end(); ++it) {
            if (((*it)->getPolygon())->isPntInPolygon (*station)) {
                bool rval(dynamic_cast<GMSHPolygonTree*>((*it))->insertStation (station));
                // stop recursion if sub SimplePolygonTree is a leaf
                if (rval && (*it)->getNumberOfChildren() == 0)
                    _stations.push_back (station);
                return rval;
            }
        }
        // station did not fit into child nodes -> insert the station into this node
        _stations.push_back (station);
        return true;
    } else {
        return false;
    }
}

void GMSHPolygonTree::insertPolyline(GeoLib::PolylineWithSegmentMarker * ply)
{
    if (!_node_polygon->isPartOfPolylineInPolygon(*ply))
        return;

    // check if polyline segments are inside of the polygon, intersect the
    // polygon or are part of the boundary of the polygon
    for (auto * polygon_tree : _children) {
        dynamic_cast<GMSHPolygonTree*>(polygon_tree)->insertPolyline(ply);
    }

    // calculate possible intersection points between the node polygon
    // (_node_polygon) and the given polyline ply
    // pay attention: loop bound is not fix!
    GeoLib::Point tmp_pnt;
    GeoLib::PointVec & pnt_vec(*(_geo_objs.getPointVecObj(_geo_name)));
    for (auto segment_it(ply->begin()); segment_it != ply->end();
         ++segment_it)
    {
        if (ply->isSegmentMarked(segment_it.getSegmentNumber()))
            continue;

        if (_node_polygon->containsSegment(*segment_it)) {
            ply->markSegment(segment_it.getSegmentNumber(), true);
            continue;
        }

        std::size_t seg_num(0);
        GeoLib::Point intersection_pnt;
        while (_node_polygon->getNextIntersectionPointPolygonLine(
            *segment_it, intersection_pnt, seg_num))
        {
            // insert the intersection point to point vector of GEOObjects instance
            const std::size_t pnt_vec_size(pnt_vec.size());
            std::size_t pnt_id(
                pnt_vec.push_back(new GeoLib::Point(intersection_pnt)));
            if (pnt_vec_size < pnt_vec.size()) { // case: new point
                // modify the polygon
                _node_polygon->insertPoint(seg_num+1, pnt_id);
                // modify the polyline
                ply->insertPoint(segment_it.getSegmentNumber(), pnt_id);
            } else { // case: existing point
                // check if point id is within the polygon
                if (! _node_polygon->isPointIDInPolyline(pnt_id)) {
                    _node_polygon->insertPoint(seg_num+1, pnt_id);
                }

                // check if point id is in polyline
                if (! ply->isPointIDInPolyline(pnt_id)) {
                    ply->insertPoint(segment_it.getSegmentNumber()+1, pnt_id);
                }
            }

            std::size_t tmp_seg_num(seg_num+1);
            if (!_node_polygon->getNextIntersectionPointPolygonLine(
                    *segment_it, tmp_pnt, tmp_seg_num))
            {
                // check a point of the segment except the end points
                for (std::size_t i(0); i<3; i++) {
                    tmp_pnt[i] = ((*segment_it).getBeginPoint()[i] +
                                  (*segment_it).getEndPoint()[i]) /
                                 2;
                }
                if (_node_polygon->isPntInPolygon(tmp_pnt)) {
                    ply->markSegment(segment_it.getSegmentNumber(), true);
                    // insert line segment as constraint
                    _gmsh_lines_for_constraints.push_back(
                        new GMSHLine((*segment_it).getBeginPoint().getID(),
                                     (*segment_it).getEndPoint().getID()));
                }
            }
            seg_num++;

            // check a point of the segment except the end points
            for (std::size_t i(0); i<3; i++) {
                tmp_pnt[i] = ((*segment_it).getBeginPoint()[i] +
                              (*segment_it).getEndPoint()[i]) /
                             2;
            }

            checkIntersectionsSegmentExistingPolylines(ply, segment_it);

            if (_node_polygon->isPntInPolygon(tmp_pnt)) {
                ply->markSegment(segment_it.getSegmentNumber(), true);
                // insert line segment as constraint
                _gmsh_lines_for_constraints.push_back(
                    new GMSHLine((*segment_it).getBeginPoint().getID(),
                                 (*segment_it).getEndPoint().getID()));
            }
        }
    }

    _plys.push_back(ply);
}

void GMSHPolygonTree::checkIntersectionsSegmentExistingPolylines(
    GeoLib::PolylineWithSegmentMarker* ply,
    GeoLib::Polyline::SegmentIterator const& seg_it)
{
    std::size_t const ply_segment_number(seg_it.getSegmentNumber());
    for(GeoLib::PolylineWithSegmentMarker *const p : _plys) {
        std::size_t n_segments(p->getNumberOfSegments());
        GeoLib::PointVec & pnt_vec(*(_geo_objs.getPointVecObj(_geo_name)));
        for (auto seg_it_p(p->begin()); seg_it_p != p->end(); ++seg_it_p) {
            GeoLib::Point s; // intersection point
            if (GeoLib::lineSegmentIntersect(*seg_it, *seg_it_p, s))
            {
                const std::size_t pnt_vec_size(pnt_vec.size());
                // point id of new point in GEOObjects instance
                const std::size_t pnt_id(pnt_vec.push_back(new GeoLib::Point(s)));
                if (pnt_vec_size < pnt_vec.size()) { // case: new point
                    // modify polyline already in this node
                    p->insertPoint(seg_it_p.getSegmentNumber()+1, pnt_id);
                    n_segments++;
                    // modify polyline
                    ply->insertPoint(ply_segment_number+1, pnt_id);
                } else { // case: point exists already in geometry
                    // check if point is not alread in polyline p
                    std::size_t const k(seg_it_p.getSegmentNumber());
                    if (p->getPointID(k) != pnt_id && p->getPointID(k+1) != pnt_id) {
                        p->insertPoint(k+1, pnt_id);
                        n_segments++;
                    }
                    // check if point is not already in polyline ply
                    if (ply->getPointID(ply_segment_number) != pnt_id
                            && ply->getPointID(ply_segment_number+1) != pnt_id) {
                        ply->insertPoint(ply_segment_number+1, pnt_id);
                    }
                }
            }
        }
    }
}

void GMSHPolygonTree::initMeshDensityStrategy()
{
    if (dynamic_cast<GMSHAdaptiveMeshDensity*> (_mesh_density_strategy)) {
        // collect points
        std::vector<GeoLib::Point const*> pnts;
        const std::size_t n_pnts_polygon (_node_polygon->getNumberOfPoints());
        for (std::size_t k(0); k<n_pnts_polygon; k++) {
            pnts.push_back(_node_polygon->getPoint(k));
        }
        getPointsFromSubPolygons(pnts);

        const std::size_t n_plys (_plys.size());
        for (std::size_t k(0); k<n_plys; k++) {
            const std::size_t n_pnts_in_kth_ply(_plys[k]->getNumberOfPoints());
            for (std::size_t j(0); j<n_pnts_in_kth_ply; j++) {
                pnts.push_back(_plys[k]->getPoint(j));
            }
        }

        // give collected points to the mesh density strategy
        _mesh_density_strategy->initialize(pnts);
        // insert constraints
        dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)->addPoints(_stations);
        std::vector<GeoLib::Point const*> stations;
        getStationsInsideSubPolygons(stations);
        dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)->addPoints(stations);
    }
}

void GMSHPolygonTree::createGMSHPoints(std::vector<GMSHPoint*> & gmsh_pnts) const
{
    const std::size_t n_pnts_polygon (_node_polygon->getNumberOfPoints());
    for (std::size_t k(0); k<n_pnts_polygon-1; k++) {
        const std::size_t id (_node_polygon->getPointID(k));
        GeoLib::Point const*const pnt(_node_polygon->getPoint(k));
        // if this point was already part of another polyline
        if (gmsh_pnts[id] != nullptr)
            continue;
        gmsh_pnts[id] = new GMSHPoint(
            *pnt, id, _mesh_density_strategy->getMeshDensityAtPoint(pnt));
    }

    const std::size_t n_plys(_plys.size());
    for (std::size_t k(0); k<n_plys; k++) {
        const std::size_t n_pnts_in_ply(_plys[k]->getNumberOfPoints());
        for (std::size_t j(0); j<n_pnts_in_ply; j++) {
            if (_node_polygon->isPntInPolygon(*(_plys[k]->getPoint(j)))) {
                const std::size_t id (_plys[k]->getPointID(j));
                // if this point was already part of another polyline
                if (gmsh_pnts[id] != nullptr)
                    continue;
                GeoLib::Point const*const pnt(_plys[k]->getPoint(j));
                gmsh_pnts[id] = new GMSHPoint(
                    *pnt, id,
                    _mesh_density_strategy->getMeshDensityAtPoint(pnt));
            }
        }
    }

    // walk through children
    for (auto it : _children)
    {
        dynamic_cast<GMSHPolygonTree*>(it)->createGMSHPoints(gmsh_pnts);
    }
}

void GMSHPolygonTree::writeLineLoop(std::size_t &line_offset, std::size_t &sfc_offset, std::ostream& out) const
{
    const std::size_t n_pnts (_node_polygon->getNumberOfPoints());
    std::size_t first_pnt_id(_node_polygon->getPointID(0)), second_pnt_id;
    for (std::size_t k(1); k<n_pnts; k++) {
        second_pnt_id = _node_polygon->getPointID(k);
        out << "Line(" << line_offset + k-1 << ") = {" << first_pnt_id << "," << second_pnt_id << "};\n";
        first_pnt_id = second_pnt_id;
    }
    out << "Line Loop(" << line_offset + n_pnts-1 << ") = {";
    for (std::size_t k(0); k<n_pnts - 2; k++) {
        out << line_offset+k << ",";
    }
    out << line_offset+n_pnts-2 << "};\n";
    out << "Plane Surface(" << sfc_offset << ") = {" << line_offset+n_pnts-1 << "};\n";
    line_offset += n_pnts;
    sfc_offset++;
}

void GMSHPolygonTree::writeLineConstraints(std::size_t& line_offset,
                                           std::size_t sfc_number,
                                           std::ostream& out) const
{
    for (auto polyline : _plys)
    {
        const std::size_t n_pnts(polyline->getNumberOfPoints());
        std::size_t first_pnt_id(polyline->getPointID(0)), second_pnt_id;
        for (std::size_t k(1); k < n_pnts; k++)
        {
            second_pnt_id = polyline->getPointID(k);
            if (polyline->isSegmentMarked(k - 1) &&
                _node_polygon->isPntInPolygon(*(polyline->getPoint(k))) &&
                !GeoLib::containsEdge(*_node_polygon, first_pnt_id, second_pnt_id))
            {
                out << "Line(" << line_offset + k - 1 << ") = {" << first_pnt_id
                    << "," << second_pnt_id << "};\n";
                out << "Line { " << line_offset + k - 1 << " } In Surface { "
                    << sfc_number << " };\n";
            }
            first_pnt_id = second_pnt_id;
        }
        line_offset += n_pnts;
    }
}

void GMSHPolygonTree::writeSubPolygonsAsLineConstraints(std::size_t &line_offset, std::size_t sfc_number, std::ostream& out) const
{
    for (auto it : _children)
    {
        dynamic_cast<GMSHPolygonTree*>(it)->writeSubPolygonsAsLineConstraints(
            line_offset, sfc_number, out);
    }

    if (_parent != nullptr)
    {
        const std::size_t n_pnts(_node_polygon->getNumberOfPoints());
        std::size_t first_pnt_id(_node_polygon->getPointID(0)), second_pnt_id;
        for (std::size_t k(1); k<n_pnts; k++) {
            second_pnt_id = _node_polygon->getPointID(k);
            out << "Line(" << line_offset + k-1 << ") = {" << first_pnt_id << "," << second_pnt_id << "};\n";
            first_pnt_id = second_pnt_id;
            out << "Line { " << line_offset+k-1 << " } In Surface { " << sfc_number << " };\n";
        }
        line_offset += n_pnts;
    }

}

void GMSHPolygonTree::writeStations(std::size_t & pnt_id_offset, std::size_t sfc_number, std::ostream& out) const
{
    for (auto const* station : _stations) {
        out << "Point(" << pnt_id_offset << ") = {" << (*station)[0] << ", "
            << (*station)[1] << ", 0.0, "
            << _mesh_density_strategy->getMeshDensityAtStation(station)
            << "}; // Station "
            << static_cast<GeoLib::Station const*>(station)->getName() << " \n";
        out << "Point { " << pnt_id_offset << " } In Surface { " << sfc_number << " };\n";
        ++pnt_id_offset;
    }
}

void GMSHPolygonTree::writeAdditionalPointData(std::size_t & pnt_id_offset, std::size_t sfc_number, std::ostream& out) const
{
    if (dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)) {
        std::vector<GeoLib::Point*> steiner_pnts;
        dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)->getSteinerPoints(steiner_pnts, 0);
        const std::size_t n(steiner_pnts.size());
        for (std::size_t k(0); k<n; k++) {
            if (_node_polygon->isPntInPolygon(*(steiner_pnts[k]))) {
                out << "Point(" << pnt_id_offset + k << ") = {" << (*(steiner_pnts[k]))[0] << "," << (*(steiner_pnts[k]))[1] << ", 0.0, ";
                out << _mesh_density_strategy->getMeshDensityAtPoint(steiner_pnts[k]) << "};\n";
                out << "Point { " << pnt_id_offset + k << " } In Surface { " << sfc_number << " };\n";
            }
            delete steiner_pnts[k];
        }
        pnt_id_offset += n;
    }

#ifndef NDEBUG
    if (dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)) {
        auto pnts = std::unique_ptr<std::vector<GeoLib::Point*>>(
            new std::vector<GeoLib::Point*>);
        auto plys = std::unique_ptr<std::vector<GeoLib::Polyline*>>(
            new std::vector<GeoLib::Polyline*>);
        dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)->getQuadTreeGeometry(*pnts, *plys);
        std::string quad_tree_geo("QuadTree");
        _geo_objs.addPointVec(std::move(pnts), quad_tree_geo);
        std::vector<std::size_t> const& id_map ((_geo_objs.getPointVecObj(quad_tree_geo))->getIDMap());
        for (std::size_t k(0); k<plys->size(); k++) {
            for (std::size_t j(0); j<(*plys)[k]->getNumberOfPoints(); j++) {
                ((*plys)[k])->setPointID(j, id_map[((*plys)[k])->getPointID(j)]);
            }
        }
        _geo_objs.addPolylineVec(std::move(plys), quad_tree_geo);
    }
#endif

}

void GMSHPolygonTree::getPointsFromSubPolygons(std::vector<GeoLib::Point const*>& pnts)
{
    for (std::list<SimplePolygonTree*>::const_iterator it (_children.begin()); it != _children.end(); ++it) {
        dynamic_cast<GMSHPolygonTree*>((*it))->getPointsFromSubPolygons(pnts);
    }
}

void GMSHPolygonTree::getStationsInsideSubPolygons(std::vector<GeoLib::Point const*>& stations)
{
    const std::size_t n_stations(_stations.size());
    for (std::size_t k(0); k<n_stations; k++) {
        stations.push_back(_stations[k]);
    }

    for (std::list<SimplePolygonTree*>::const_iterator it (_children.begin()); it != _children.end(); ++it) {
        dynamic_cast<GMSHPolygonTree*>((*it))->getStationsInsideSubPolygons(stations);
    }
}

}  // end namespace GMSH
}  // end namespace FileIO
