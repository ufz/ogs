/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <memory>

#include "gtest/gtest.h"

#include "BaseLib/BuildInfo.h"
#include "Applications/FileIO/SWMM/SwmmInterface.h"
#include "GeoLib/GeoObjects.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshInformation.h"

TEST(FileIO, TestSwmmInterface)
{
    std::string const base_name ("swmm_input_example");
    std::string const file_name (BaseLib::BuildInfo::data_path + "/FileIO/" + base_name + ".inp");

    // testing geo conversion
    GeoLib::GEOObjects geo_objects;
    bool const result =
        FileIO::SwmmInterface::convertSwmmInputToGeometry(file_name, geo_objects, true);
    ASSERT_TRUE(result);

    std::vector<GeoLib::Point*> const*const pnts (geo_objects.getPointVec(base_name));
    std::size_t const n_pnts (259);
    ASSERT_EQ(n_pnts, pnts->size());

    for (std::size_t i=0; i<21; i++)
        ASSERT_TRUE((*(*pnts)[i])[2] > 0);

    for (std::size_t i=21; i<n_pnts; i++)
        ASSERT_NEAR(0, (*(*pnts)[i])[2], std::numeric_limits<double>::epsilon());

    std::vector<GeoLib::Polyline*> const*const lines (geo_objects.getPolylineVec(base_name));
    std::size_t n_lines (79);
    ASSERT_EQ(n_lines, lines->size());

    for (std::size_t i=0; i<58; ++i)
        ASSERT_TRUE((*lines)[i]->getNumberOfPoints() > 3 && (*lines)[i]->getNumberOfPoints() < 7);
    for (std::size_t i=58; i<n_lines; ++i)
        ASSERT_TRUE(2 == (*lines)[i]->getNumberOfPoints());

    // testing mesh conversion
    std::unique_ptr<FileIO::SwmmInterface> swmm = FileIO::SwmmInterface::create(file_name);
    ASSERT_TRUE(swmm != nullptr);

    MeshLib::Mesh const& mesh (swmm->getMesh());
    std::size_t const n_nodes (22);
    std::size_t const n_elems (21);
    ASSERT_EQ(n_nodes, mesh.getNumberOfNodes());
    ASSERT_EQ(n_elems, mesh.getNumberOfElements());

    auto const* const depth =
        mesh.getProperties().getPropertyVector<double>("Max Depth");
    ASSERT_TRUE(n_nodes == depth->size());
    ASSERT_NEAR((*depth)[1], 2*(*depth)[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR((*depth)[2], 2*(*depth)[1], std::numeric_limits<double>::epsilon());
    for (std::size_t i=3; i<n_nodes; ++i)
        ASSERT_NEAR(0, (*depth)[i], std::numeric_limits<double>::epsilon());

    std::array<unsigned, 7> types (MeshLib::MeshInformation::getNumberOfElementTypes(mesh));
    ASSERT_EQ(n_elems, types[0]); // all elems are lines
    std::pair<int, int> bounds (MeshLib::MeshInformation::getValueBounds<int>(mesh, "MaterialIDs"));
    ASSERT_EQ(0, bounds.first);
    ASSERT_EQ(0, bounds.second);

    ASSERT_NEAR(186.06,  mesh.getMinEdgeLength(), 0.01);
    ASSERT_NEAR(3171.42, mesh.getMaxEdgeLength(), 0.01);
}

