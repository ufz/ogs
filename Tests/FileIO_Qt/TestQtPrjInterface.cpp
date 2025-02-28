/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <gtest/gtest.h>

#include <cstdio>

#include "Applications/DataHolderLib/BoundaryCondition.h"
#include "Applications/DataHolderLib/Project.h"
#include "Applications/DataHolderLib/SourceTerm.h"
#include "Applications/FileIO/XmlIO/Qt/XmlPrjInterface.h"
#include "GeoLib/GEOObjects.h"
#include "InfoLib/TestInfo.h"

TEST(TestQtPrjInterface, QtXmlPrjReader)
{
    struct TestParams
    {
        std::string const file_name;
        std::size_t n_geo;
        std::size_t n_mesh;
        std::size_t n_bc;
        std::size_t n_st;
    };

    std::vector<TestParams> test_files;

    std::string name =
        TestInfoLib::TestInfo::data_path +
        "/Elliptic/nonuniform_bc_SteadyStateDiffusion/neumann_nonuniform.prj";
    test_files.push_back({name, 0, 3, 2, 0});
    name =
        TestInfoLib::TestInfo::data_path +
        "/Elliptic/nonuniform_bc_SteadyStateDiffusion/neumann_nonuniform.prj";
    test_files.push_back({name, 0, 3, 2, 0});

    for (auto& test_file : test_files)
    {
        DataHolderLib::Project project;
        FileIO::XmlPrjInterface xml(project);
        int result = xml.readFile(test_file.file_name);
        EXPECT_EQ(1, result);

        auto const geo_names = project.getGEOObjects().getGeometryNames();
        EXPECT_EQ(test_file.n_geo, geo_names.size());
        EXPECT_EQ(test_file.n_mesh, project.getMeshObjects().size());

        std::vector<std::unique_ptr<DataHolderLib::BoundaryCondition>> const&
            bcs = project.getBoundaryConditions();
        EXPECT_EQ(test_file.n_bc, bcs.size());
        for (auto& bc : bcs)
        {
            EXPECT_FALSE(bc->getProcessVarName().empty());
            EXPECT_FALSE(bc->getParamName().empty());
            EXPECT_FALSE(
                DataHolderLib::BoundaryCondition::ConditionType::NONE ==
                bc->getType());
            EXPECT_FALSE(bc->getBaseObjName().empty());
            if (bc->getBaseObjType() == DataHolderLib::BaseObjType::GEOMETRY)
            {
                EXPECT_FALSE(bc->getObjName().empty());
            }
        }

        std::vector<std::unique_ptr<DataHolderLib::SourceTerm>> const& sts =
            project.getSourceTerms();
        EXPECT_EQ(test_file.n_st, sts.size());
        for (auto& st : sts)
        {
            EXPECT_FALSE(st->getProcessVarName().empty());
            EXPECT_FALSE(st->getParamName().empty());
            EXPECT_FALSE(DataHolderLib::SourceTerm::ConditionType::NONE ==
                         st->getType());
            EXPECT_FALSE(st->getBaseObjName().empty());
            if (st->getBaseObjType() == DataHolderLib::BaseObjType::GEOMETRY)
            {
                EXPECT_FALSE(st->getObjName().empty());
            }
        }
    }
}
