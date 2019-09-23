/**
 * @brief Tests for GeoLib::GEOOnjects::geoPointsToStations()
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "InfoLib/TestInfo.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"

TEST(GeoLib, PointToStationConversion)
{
    std::string const file_name(
        TestInfoLib::TestInfo::data_path +
        "/MeshGeoToolsLib/Hamburg/BrunnenGOK_Outline_ModellHH.gml");
    GeoLib::GEOObjects geo_obj;
    GeoLib::IO::BoostXmlGmlInterface io(geo_obj);
    bool const ret = io.readFile(file_name);
    EXPECT_TRUE(ret);
    std::vector<std::string> geo_names;
    geo_obj.getGeometryNames(geo_names);
    auto const* const pnts = geo_obj.getPointVec(geo_names[0]);
    assert(pnts != nullptr);
    std::size_t const exp_all_pnts(310);
    EXPECT_EQ(exp_all_pnts, pnts->size());

    // converting only unused points
    std::string stn_orphaned_pnts("Orphaned Points");
    int const res_orphaned_pnts =
        geo_obj.geoPointsToStations(geo_names[0], stn_orphaned_pnts, true);
    EXPECT_EQ(0, res_orphaned_pnts);
    auto const* const stns = geo_obj.getStationVec(stn_orphaned_pnts);
    assert(stns != nullptr);
    std::size_t const n_stns = stns->size();
    EXPECT_EQ(18, n_stns);
    for (std::size_t i = 0; i < n_stns; ++i)
    {
        EXPECT_EQ((*(*stns)[i])[0], (*(*pnts)[i])[0]);
        EXPECT_EQ((*(*stns)[i])[1], (*(*pnts)[i])[1]);
        EXPECT_EQ((*(*stns)[i])[2], (*(*pnts)[i])[2]);
    }

    // converting all points
    std::string stn_all_pnts("All Points");
    int const res_all_pnts =
        geo_obj.geoPointsToStations(geo_names[0], stn_all_pnts, false);
    EXPECT_EQ(0, res_all_pnts);
    EXPECT_EQ(exp_all_pnts, geo_obj.getStationVec(stn_all_pnts)->size());
}
