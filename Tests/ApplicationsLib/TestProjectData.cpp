/**
 * \file   TestProjectData.cpp
 * \date   2015-07-13
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include "gtest/gtest.h"
#include "Applications/ApplicationsLib/ProjectData.h"

TEST(ApplicationsLib, ProjectData)
{
    ProjectData project;
    GeoLib::GEOObjects *geo_objects = project.getGEOObjects();

    ASSERT_TRUE(geo_objects != nullptr);
}
