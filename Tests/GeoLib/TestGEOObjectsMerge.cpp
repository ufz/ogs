/**
 * @file TestGEOObjectsMerge.cpp
 * @author Thomas Fischer
 * @date May 21, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <map>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "GeoLib/GEOObjects.h"

void createSetOfTestPointsAndAssociatedNames(GeoLib::GEOObjects & geo_objs, std::string &name, GeoLib::Point const& shift)
{
    auto pnts = std::unique_ptr<std::vector<GeoLib::Point*>>(
        new std::vector<GeoLib::Point*>);
    auto* pnt_name_map(new std::map<std::string, std::size_t>);

    const std::size_t pnts_per_edge(8);
    for (std::size_t k(0); k < pnts_per_edge; k++) {
        const std::size_t k_offset(k * pnts_per_edge * pnts_per_edge);
        for (std::size_t j(0); j < pnts_per_edge; j++) {
            const std::size_t offset(j * pnts_per_edge + k_offset);
            for (std::size_t i(0); i < pnts_per_edge; i++) {
                std::size_t const id(i+offset);
                pnts->push_back(
                    new GeoLib::Point(i+shift[0], j+shift[1], k+shift[2], id));
                std::string pnt_name(
                        name + "-" + std::to_string(i) + "-" + std::to_string(j) + "-"
                                + std::to_string(k));
                pnt_name_map->insert(std::make_pair(pnt_name, id));
            }
        }
    }

    geo_objs.addPointVec(std::move(pnts), name, pnt_name_map);
}

TEST(GeoLib, GEOObjectsMergePoints)
{
    GeoLib::GEOObjects geo_objs;
    std::vector<std::string> names;

    // *** insert set of points number 0
    GeoLib::Point shift (0.0,0.0,0.0);
    names.emplace_back("PointSet0");
    createSetOfTestPointsAndAssociatedNames(geo_objs, names[0], shift);

    // *** insert set of points number 1
    names.emplace_back("PointSet1");
    createSetOfTestPointsAndAssociatedNames(geo_objs, names[1], shift);

    // *** merge geometries
    std::string merged_geometries_name("MergedEqualPointSet");
    geo_objs.mergeGeometries(names, merged_geometries_name);

    GeoLib::PointVec const* merged_point_vec (geo_objs.getPointVecObj(merged_geometries_name));

    ASSERT_TRUE(merged_point_vec != nullptr);
    ASSERT_EQ(512u, merged_point_vec->size());
    std::string test_name;
    merged_point_vec->getNameOfElementByID(0, test_name);
    ASSERT_EQ("PointSet0-0-0-0", test_name);
    merged_point_vec->getNameOfElementByID(511, test_name);
    ASSERT_EQ("PointSet0-7-7-7", test_name);

    // *** insert "shifted" set of points
    shift[0] += 1e-4;
    names.emplace_back("ShiftedPointSet");
    createSetOfTestPointsAndAssociatedNames(geo_objs, names[2], shift);

    // *** merge PointSet0, PointSet1 and ShiftedPointSet
    merged_geometries_name = "MergedShiftedPointSet";
    geo_objs.mergeGeometries(names, merged_geometries_name);
    merged_point_vec = geo_objs.getPointVecObj(merged_geometries_name);

    ASSERT_TRUE(merged_point_vec != nullptr);
    ASSERT_EQ(1024u, merged_point_vec->size());
    merged_point_vec->getNameOfElementByID(0, test_name);
    ASSERT_EQ("PointSet0-0-0-0", test_name);
    merged_point_vec->getNameOfElementByID(511, test_name);
    ASSERT_EQ("PointSet0-7-7-7", test_name);
    merged_point_vec->getNameOfElementByID(512, test_name);
    ASSERT_EQ("ShiftedPointSet-0-0-0", test_name);
    merged_point_vec->getNameOfElementByID(1023, test_name);
    ASSERT_EQ("ShiftedPointSet-7-7-7", test_name);

    std::size_t id;
    ASSERT_TRUE(merged_point_vec->getElementIDByName (test_name, id));
    ASSERT_EQ(1023u, id);

    test_name = "PointSet1-0-0-0";
    ASSERT_FALSE(merged_point_vec->getElementIDByName (test_name, id));
}

TEST(GeoLib, GEOObjectsMergePointsAndPolylines)
{
    GeoLib::GEOObjects geo_objs;
    std::vector<std::string> names;

    // *** insert points to vector
    auto pnts = std::unique_ptr<std::vector<GeoLib::Point*>>(
        new std::vector<GeoLib::Point*>);
    pnts->reserve(4);
    pnts->push_back(new GeoLib::Point(0.0,0.0,0.0));
    pnts->push_back(new GeoLib::Point(1.0,0.0,0.0));
    pnts->push_back(new GeoLib::Point(1.0,1.0,0.0));
    pnts->push_back(new GeoLib::Point(0.0,1.0,0.0));

    std::string geometry_0("GeometryWithPntsAndPolyline");
    geo_objs.addPointVec(std::move(pnts), geometry_0, nullptr, std::numeric_limits<double>::epsilon());

    // *** insert polyline
    auto* ply(new GeoLib::Polyline(*geo_objs.getPointVec(geometry_0)));
    ply->addPoint(0);
    ply->addPoint(1);
    ply->addPoint(2);
    ply->addPoint(3);
    ply->addPoint(0);
    auto plys = std::unique_ptr<std::vector<GeoLib::Polyline*>>(
        new std::vector<GeoLib::Polyline*>);
    plys->push_back(ply);
    geo_objs.addPolylineVec(std::move(plys), geometry_0, nullptr);
    names.push_back(geometry_0);

    // *** insert set of points number
    GeoLib::Point shift (0.0,0.0,0.0);
    names.emplace_back("PointSet0");
    createSetOfTestPointsAndAssociatedNames(geo_objs, names[1], shift);

    // *** merge geometries
    std::string merged_geometries_name("MergedQuadGeoAndPointSet");
    geo_objs.mergeGeometries(names, merged_geometries_name);

    std::vector<GeoLib::Polyline*> const* const polylines =
        geo_objs.getPolylineVec(merged_geometries_name);

    ASSERT_TRUE(polylines != nullptr);
    ASSERT_EQ(1u, polylines->size());
}

TEST(GeoLib, GEOObjectsMergePolylinesWithNames)
{
    GeoLib::GEOObjects geo_objs;
    std::vector<std::string> names;

    // *** insert first set of points to vector (for first polyline)
    auto pnts_0 = std::unique_ptr<std::vector<GeoLib::Point*>>(
        new std::vector<GeoLib::Point*>);
    pnts_0->reserve(4);
    pnts_0->push_back(new GeoLib::Point(0.0,0.0,0.0));
    pnts_0->push_back(new GeoLib::Point(1.0,0.0,0.0));
    pnts_0->push_back(new GeoLib::Point(1.0,1.0,0.0));
    pnts_0->push_back(new GeoLib::Point(0.0,1.0,0.0));

    std::string geometry_0("Geometry0");
    geo_objs.addPointVec(std::move(pnts_0),
                         geometry_0,
                         nullptr,
                         std::numeric_limits<double>::epsilon());

    // *** insert a named polyline into geometry
    auto* ply_00(new GeoLib::Polyline(*geo_objs.getPointVec(geometry_0)));
    ply_00->addPoint(0);
    ply_00->addPoint(1);
    ply_00->addPoint(2);
    ply_00->addPoint(3);
    ply_00->addPoint(0);
    auto plys_0 = std::unique_ptr<std::vector<GeoLib::Polyline*>>(
        new std::vector<GeoLib::Polyline*>);
    plys_0->push_back(ply_00);
    auto* names_map_0(new std::map<std::string, std::size_t>);
    names_map_0->insert(std::pair<std::string, std::size_t>("Polyline0FromGeometry0", 0));
    geo_objs.addPolylineVec(std::move(plys_0), geometry_0, names_map_0);
    names.push_back(geometry_0);

    auto pnts_1 = std::unique_ptr<std::vector<GeoLib::Point*>>(
        new std::vector<GeoLib::Point*>);
    pnts_1->reserve(4);
    pnts_1->push_back(new GeoLib::Point(0.0,0.0,0.0));
    pnts_1->push_back(new GeoLib::Point(1.0,0.0,0.0));
    pnts_1->push_back(new GeoLib::Point(1.0,1.0,0.0));
    pnts_1->push_back(new GeoLib::Point(0.0,1.0,0.0));

    std::string geometry_1("Geometry1");
    geo_objs.addPointVec(std::move(pnts_1),
                         geometry_1,
                         nullptr,
                         std::numeric_limits<double>::epsilon());

    // *** insert a named polyline into geometry
    auto* ply_10(new GeoLib::Polyline(*geo_objs.getPointVec(geometry_1)));
    ply_10->addPoint(0);
    ply_10->addPoint(1);
    auto* ply_11(new GeoLib::Polyline(*geo_objs.getPointVec(geometry_1)));
    ply_11->addPoint(2);
    ply_11->addPoint(3);
    auto plys_1 = std::unique_ptr<std::vector<GeoLib::Polyline*>>(
        new std::vector<GeoLib::Polyline*>);
    plys_1->push_back(ply_10);
    plys_1->push_back(ply_11);
    auto* names_map_1(new std::map<std::string, std::size_t>);
    names_map_1->insert(std::pair<std::string, std::size_t>("Polyline0FromGeometry1", 0));
    names_map_1->insert(std::pair<std::string, std::size_t>("Polyline1FromGeometry1", 1));
    geo_objs.addPolylineVec(std::move(plys_1), geometry_1, names_map_1);
    names.push_back(geometry_1);

    // *** merge geometries
    std::string merged_geometries_name("MergedPolylinesWithNames");
    geo_objs.mergeGeometries(names, merged_geometries_name);

    // *** tests
    // check number of points
    ASSERT_EQ(4u, geo_objs.getPointVec(merged_geometries_name)->size());

    GeoLib::PolylineVec const*const ply_vec_objs =
        geo_objs.getPolylineVecObj(merged_geometries_name);
    std::vector<GeoLib::Polyline*> const* const polylines =
        ply_vec_objs->getVector();

    // check number of polylines
    ASSERT_TRUE(polylines != nullptr);
    ASSERT_EQ(3u, polylines->size());

    // check names of polylines
    ASSERT_TRUE(ply_vec_objs->getElementByName("Polyline0FromGeometry0") != nullptr);
    ASSERT_TRUE(ply_vec_objs->getElementByName("Polyline0FromGeometry1") != nullptr);
    ASSERT_TRUE(ply_vec_objs->getElementByName("Polyline1FromGeometry1") != nullptr);
}

