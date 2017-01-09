/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <cstdlib>
#include <ctime>
#include <random>
#include <string>

#include "gtest/gtest.h"

#include "GeoLib/DuplicateGeometry.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"
#include "GeoLib/Triangle.h"

TEST(GeoLib, DuplicateGeometry)
{
    GeoLib::GEOObjects geo;
    std::string input_name ("input");

    // generate points
    std::srand ( static_cast<unsigned>(std::time(nullptr)) );
    std::size_t n_pnts (rand() % 1000 + 100);
    int box_size (std::rand());
    double half_box_size(box_size/2);
    std::unique_ptr<std::vector<GeoLib::Point*>> pnts (new std::vector<GeoLib::Point*>);
    pnts->reserve(n_pnts);
    for (int k(0); k<static_cast<int>(n_pnts); k++) {
        pnts->push_back(new GeoLib::Point(
            std::rand() % box_size - half_box_size,
            std::rand() % box_size - half_box_size,
            std::rand() % box_size - half_box_size));
    }
    geo.addPointVec(std::move(pnts), input_name);
    // update number of points in case because possibly non-unique points have been removed
    n_pnts = geo.getPointVec(input_name)->size();
    std::string output ("output_geometry");

    // duplicate points
    {
        GeoLib::DuplicateGeometry dup (geo, input_name, output);

        std::vector<GeoLib::Point*> const*const pnts (geo.getPointVec(input_name));
        std::vector<GeoLib::Point*> const*const new_pnts (geo.getPointVec(output));
        std::vector<GeoLib::Point*>& mod_pnts (dup.getPointVectorCopy());
        ASSERT_EQ(n_pnts, new_pnts->size());
        ASSERT_EQ(n_pnts, mod_pnts.size());

        for (std::size_t i=0; i<n_pnts; ++i)
        {
            ASSERT_EQ((*(*pnts)[i])[0], (*(*new_pnts)[i])[0]);
            ASSERT_EQ((*(*pnts)[i])[1], (*(*new_pnts)[i])[1]);
            ASSERT_EQ((*(*pnts)[i])[2], (*(*new_pnts)[i])[2]);
        }
        mod_pnts.push_back(new GeoLib::Point(0,0,0));
        ASSERT_EQ(mod_pnts.size(), pnts->size() + 1);
        ASSERT_EQ(mod_pnts.size(), new_pnts->size());
    }

    // add polylines
    std::size_t n_plys (rand() % 10 + 1);
    std::unique_ptr<std::vector<GeoLib::Polyline*>> plys (new std::vector<GeoLib::Polyline*>);
    for (std::size_t i=0; i<n_plys; ++i)
    {
        int n_ply_pnts (rand() % 100 + 2);
        GeoLib::Polyline* line = new GeoLib::Polyline(*geo.getPointVec(input_name));
        for (std::size_t j=0; j<static_cast<std::size_t>(n_ply_pnts); ++j)
            line->addPoint(rand() % n_pnts);
        plys->push_back(line);
    }
    geo.addPolylineVec(std::move(plys), input_name);

    // duplicate polylines
    {
        GeoLib::DuplicateGeometry dup (geo, input_name, output);
        std::string const& output2 = dup.getFinalizedOutputName();
        ASSERT_FALSE(output == output2);

        std::vector<GeoLib::Polyline*> const*const plys (geo.getPolylineVec(input_name));
        std::vector<GeoLib::Polyline*> const*const new_plys (geo.getPolylineVec(output2));
        ASSERT_EQ(n_plys, new_plys->size());

        for (std::size_t i=0; i<n_plys; ++i)
        {
            std::size_t const n_ply_pnts ((*new_plys)[i]->getNumberOfPoints());
            ASSERT_EQ(n_ply_pnts, (*plys)[i]->getNumberOfPoints());
            ASSERT_EQ((*new_plys)[i]->getNumberOfSegments(), (*plys)[i]->getNumberOfSegments());
            for (std::size_t j=0; j<n_ply_pnts; ++j)
                ASSERT_EQ((*plys)[i]->getPointID(j), (*new_plys)[i]->getPointID(j));
        }

        std::vector<GeoLib::Polyline*>& mod_plys (dup.getPolylineVectorCopy());
        mod_plys.push_back(new GeoLib::Polyline(*(*new_plys)[0]));
        ASSERT_EQ(mod_plys.size(), plys->size() + 1);
        ASSERT_EQ(mod_plys.size(), new_plys->size());
    }

    // add surfaces
    std::size_t n_sfcs (rand() % 10 + 1);
    std::unique_ptr<std::vector<GeoLib::Surface*>> sfcs (new std::vector<GeoLib::Surface*>);
    for (std::size_t i=0; i<n_sfcs; ++i)
    {
        int n_tris (rand() % 10);
        GeoLib::Surface* sfc = new GeoLib::Surface(*geo.getPointVec(input_name));
        for (std::size_t j=0; j<static_cast<std::size_t>(n_tris); ++j)
            sfc->addTriangle(rand() % n_pnts, rand() % n_pnts, rand() % n_pnts);
        if (sfc->getNumberOfTriangles() > 0)
            sfcs->push_back(sfc);
        else
            delete sfc;
    }
    n_sfcs = sfcs->size();
    geo.addSurfaceVec(std::move(sfcs), input_name);

    // duplicate surfaces
    {
        GeoLib::DuplicateGeometry dup (geo, input_name, output);
        std::string const& output2 = dup.getFinalizedOutputName();
        ASSERT_FALSE(output == output2);

        std::vector<GeoLib::Surface*> const* sfcs (geo.getSurfaceVec(input_name));
        std::vector<GeoLib::Surface*> const* new_sfcs (geo.getSurfaceVec(output2));
        ASSERT_EQ(n_sfcs, new_sfcs->size());

        for (std::size_t i=0; i<n_sfcs; ++i)
        {
            std::size_t const n_tris ((*new_sfcs)[i]->getNumberOfTriangles());
            ASSERT_EQ(n_tris, (*sfcs)[i]->getNumberOfTriangles());
            for (std::size_t j=0; j<n_tris; ++j)
                for (std::size_t k=0; k<3; ++k)
                    ASSERT_EQ((*(*(*sfcs)[i])[j])[k], (*(*(*new_sfcs)[i])[j])[k]);
        }

        std::vector<GeoLib::Point*>& mod_pnts (dup.getPointVectorCopy());
        std::vector<GeoLib::Surface*>& mod_sfcs (dup.getSurfaceVectorCopy());
        std::size_t n_pnts = mod_pnts.size();
        mod_pnts.push_back(new GeoLib::Point(1,0,0,n_pnts));
        mod_pnts.push_back(new GeoLib::Point(0,1,0,n_pnts+1));
        mod_pnts.push_back(new GeoLib::Point(0,0,1,n_pnts+2));
        GeoLib::Surface* sfc = new GeoLib::Surface(mod_pnts);
        sfc->addTriangle(n_pnts, n_pnts+1, n_pnts+2);
        mod_sfcs.push_back(sfc);
        ASSERT_EQ(mod_sfcs.size(), sfcs->size() + 1);
        ASSERT_EQ(mod_sfcs.size(), new_sfcs->size());
    }
}
