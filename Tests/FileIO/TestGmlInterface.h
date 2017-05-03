/**
 * \file   TestGmlInterface.h
 * \author Karsten Rink
 * \date   2016-02-21
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <cstdio>

#include "gtest/gtest.h"

#include "BaseLib/BuildInfo.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Triangle.h"

class TestGmlInterface : public testing::Test
{
public:
    GeoLib::GEOObjects geo_objects;
    std::string geo_name = "TestData";

    std::vector<GeoLib::Point> test_pnts;
    std::map<std::string, std::size_t> pnt_name_id_map;

    TestGmlInterface()
    {
        createPoints();
        createPolylines();
        createSurfaces();
    }

    void createPoints()
    {
        test_pnts.emplace_back(1,1,0,0);
        test_pnts.emplace_back(1,2,0,1);
        test_pnts.emplace_back(1,3,0,2);
        test_pnts.emplace_back(2,1,0,3);
        test_pnts.emplace_back(2,2,0,4);
        test_pnts.emplace_back(2,3,0,5);
        test_pnts.emplace_back(3,1,0,6);
        test_pnts.emplace_back(3,2,0,7);
        test_pnts.emplace_back(3,3,0,8);
        auto points = std::unique_ptr<std::vector<GeoLib::Point*>>(
            new std::vector<GeoLib::Point*>(9));

        auto cpy_name_id_map = new std::map<std::string, std::size_t>;
        std::size_t pos(0);
        for (auto p : test_pnts) {
            (*points)[pos] = new GeoLib::Point(p);
            pnt_name_id_map["p"+std::to_string(pos)] = pos;
            (*cpy_name_id_map)["p"+std::to_string(pos)] = pos;
            pos++;
        }

        geo_objects.addPointVec(std::move(points), geo_name, cpy_name_id_map);
    }

    void checkPointProperties()
    {
        auto const& pointvec(*geo_objects.getPointVecObj(geo_name));
        auto const& read_points(*pointvec.getVector());
        EXPECT_EQ(test_pnts.size(), read_points.size());
        for (std::size_t k(0); k<test_pnts.size(); ++k) {
            GeoLib::Point const& read_pnt = *read_points[k];
            // compare the coordinates
            EXPECT_EQ(test_pnts[k][0], read_pnt[0]);
            EXPECT_EQ(test_pnts[k][1], read_pnt[1]);
            EXPECT_EQ(test_pnts[k][2], read_pnt[2]);
            // compare the ids
            EXPECT_EQ(test_pnts[k].getID(), read_pnt.getID());
        }

        for (auto p : read_points) {
            // fetch name of read point
            std::string read_name;
            pointvec.getNameOfElementByID(p->getID(), read_name);
            // compare the id of the original point fetched by using the
            // read_name and the id of the read point
            EXPECT_EQ(p->getID(), pnt_name_id_map[read_name]);
        }
    }

    void createPolylines()
    {
        auto createPolyline = [](GeoLib::PointVec const& pnt_vec,
            std::vector<GeoLib::Polyline*> & lines, std::size_t line_id,
            std::vector<std::size_t> pnt_ids,
            std::map<std::string,std::size_t> & line_name_map,
            std::string const& name)
        {
            auto const & pnt_id_map(pnt_vec.getIDMap());
            lines[line_id] = new GeoLib::Polyline(*(pnt_vec.getVector()));
            for (auto pnt_id : pnt_ids)
            {
                lines[line_id]->addPoint(pnt_id_map[pnt_id]);
            }

            if (!name.empty())
                line_name_map.insert(std::make_pair(name, line_id));
        };

        auto const& pnt_vec = *(geo_objects.getPointVecObj(geo_name));

        auto lines = std::unique_ptr<std::vector<GeoLib::Polyline*>>(
            new std::vector<GeoLib::Polyline*>(5));
        auto* name_map = new std::map<std::string, std::size_t>;

        createPolyline(pnt_vec, *(lines.get()), 0, {0, 1, 2}, *name_map, "left");
        createPolyline(pnt_vec, *(lines.get()), 1, {3, 4, 5}, *name_map, "center");
        createPolyline(pnt_vec, *(lines.get()), 2, {0, 3}, *name_map, "line1");
        createPolyline(pnt_vec, *(lines.get()), 3, {3, 6}, *name_map, "line2");
        createPolyline(pnt_vec, *(lines.get()), 4, {6, 7, 8}, *name_map, "right");
        geo_objects.addPolylineVec(std::move(lines), geo_name, name_map);
    }

    void checkPolylineProperties()
    {
        const GeoLib::PolylineVec *line_vec = geo_objects.getPolylineVecObj(geo_name);
        auto const readerLines = geo_objects.getPolylineVec(geo_name);
        EXPECT_EQ(5u, readerLines->size());

        auto checkPolylineProperty = [this](GeoLib::PolylineVec const* line_vec,
            std::size_t ply_id, std::vector<std::size_t> pnt_ids,
            std::string const& name)
        {
            auto const lines = geo_objects.getPolylineVec(geo_name);
            GeoLib::Polyline* line = (*lines)[ply_id];
            EXPECT_EQ(pnt_ids.size(), line->getNumberOfPoints());
            for (std::size_t k(0); k<pnt_ids.size(); ++k)
                EXPECT_EQ(pnt_ids[k], line->getPointID(k));
            std::string line_name;
            line_vec->getNameOfElementByID(ply_id, line_name);
            EXPECT_EQ(name, line_name);
        };

        checkPolylineProperty(line_vec, 0, {0, 1, 2}, "left");
        checkPolylineProperty(line_vec, 1, {3, 4, 5}, "center");
        checkPolylineProperty(line_vec, 2, {0, 3}, "line1");
        checkPolylineProperty(line_vec, 3, {3, 6}, "line2");
        checkPolylineProperty(line_vec, 4, {6, 7, 8}, "right");
    }

    void createSurfaces()
    {
        auto const points = geo_objects.getPointVec(geo_name);
        const std::vector<std::size_t> pnt_id_map(
            geo_objects.getPointVecObj(geo_name)->getIDMap());

        auto sfcs = std::unique_ptr<std::vector<GeoLib::Surface*>>(
            new std::vector<GeoLib::Surface*>(2));
        auto* sfc_names = new std::map<std::string, std::size_t>;
        (*sfcs)[0] = new GeoLib::Surface(*points);
        (*sfcs)[0]->addTriangle(pnt_id_map[0], pnt_id_map[3], pnt_id_map[1]);
        (*sfcs)[0]->addTriangle(pnt_id_map[1], pnt_id_map[3], pnt_id_map[4]);
        (*sfcs)[0]->addTriangle(pnt_id_map[1], pnt_id_map[4], pnt_id_map[2]);
        (*sfcs)[0]->addTriangle(pnt_id_map[2], pnt_id_map[4], pnt_id_map[5]);
        (*sfc_names)["FirstSurface"] = 0;
        (*sfcs)[1] = new GeoLib::Surface(*points);
        (*sfcs)[1]->addTriangle(pnt_id_map[3], pnt_id_map[6], pnt_id_map[8]);
        (*sfcs)[1]->addTriangle(pnt_id_map[3], pnt_id_map[8], pnt_id_map[5]);
        (*sfc_names)["SecondSurface"] = 1;
        geo_objects.addSurfaceVec(std::move(sfcs), geo_name, sfc_names);
    }

    void checkSurfaceProperties()
    {
        auto checkTriangleIDs = [](
            GeoLib::Triangle const& tri, std::array<std::size_t,3> ids)
        {
            EXPECT_EQ(ids[0], tri[0]);
            EXPECT_EQ(ids[1], tri[1]);
            EXPECT_EQ(ids[2], tri[2]);
        };

        auto checkSurface = [&checkTriangleIDs](GeoLib::SurfaceVec const& sfcs,
            std::size_t sfc_id, std::vector<std::array<std::size_t,3>> tri_ids,
            std::string const& name)
        {
            auto const& sfc_vec = *(sfcs.getVector());
            auto const& sfc = *(sfc_vec[sfc_id]);
            EXPECT_EQ(tri_ids.size(), sfc.getNumberOfTriangles());
            for (std::size_t k(0); k<tri_ids.size(); ++k)
                checkTriangleIDs(*(sfc[k]), tri_ids[k]);

            std::string sfc_name;
            sfcs.getNameOfElementByID(sfc_id, sfc_name);
            EXPECT_EQ(0u, name.compare(sfc_name));
        };

        auto const read_sfcs = geo_objects.getSurfaceVecObj(geo_name);
        EXPECT_EQ(2u, read_sfcs->size());
        checkSurface(*read_sfcs, 0,
            {{{0,3,1}}, {{1,3,4}}, {{1,4,2}}, {{2,4,5}}}, "FirstSurface");
        checkSurface(*read_sfcs, 1, {{{3,6,8}}, {{3,8,5}}}, "SecondSurface");
    }
};
