/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <cstdio>

#include "gtest/gtest.h"

#include "Applications/DataHolderLib/BoundaryCondition.h"
#include "Applications/DataHolderLib/Project.h"
#include "Applications/DataHolderLib/SourceTerm.h"
#include "Applications/FileIO/XmlIO/Qt/XmlPrjInterface.h"
#include "BaseLib/BuildInfo.h"
#include "GeoLib/GEOObjects.h"

TEST(TestQtPrjInterface, QtXmlPrjReader)
{
    typedef struct
    {
        std::string const file_name;
        std::size_t n_geo;
        std::size_t n_mesh;
        std::size_t n_bc;
        std::size_t n_st;
    } TestParams;

    std::vector<TestParams> test_files;

    std::string name =
        BaseLib::BuildInfo::data_path +
        "/Elliptic/nonuniform_bc_Groundwaterflow/neumann_nonuniform.prj";
    test_files.push_back({name, 1, 1, 2, 0});
    name = BaseLib::BuildInfo::data_path +
           "/Elliptic/nonuniform_bc_Groundwaterflow/neumann_nonuniform.prj";
    test_files.push_back({name, 1, 1, 2, 0});

    for (std::size_t i = 0; i < test_files.size(); ++i)
    {
        DataHolderLib::Project project;
        FileIO::XmlPrjInterface xml(project);
        int result =
            xml.readFile(QString::fromStdString(test_files[i].file_name));
        EXPECT_EQ(1, result);

        std::vector<std::string> geo_names;
        project.getGEOObjects().getGeometryNames(geo_names);
        EXPECT_EQ(test_files[i].n_geo, geo_names.size());
        EXPECT_EQ(test_files[i].n_mesh, project.getMeshObjects().size());

        std::vector<DataHolderLib::BoundaryCondition> const bcs =
            project.getBoundaryConditions();
        EXPECT_EQ(test_files[i].n_bc, bcs.size());
        for (DataHolderLib::BoundaryCondition bc : bcs)
        {
            EXPECT_FALSE(bc.getProcessVarName().empty());
            EXPECT_FALSE(bc.getParamName().empty());
            EXPECT_FALSE(DataHolderLib::ConditionType::NONE == bc.getType());
            EXPECT_FALSE(bc.getBaseObjName().empty());
            if (bc.getBaseObjType() == DataHolderLib::BaseObjType::GEOMETRY)
                EXPECT_FALSE(bc.getObjName().empty());
        }

        std::vector<DataHolderLib::SourceTerm> const sts =
            project.getSourceTerms();
        EXPECT_EQ(test_files[i].n_st, sts.size());
        for (DataHolderLib::SourceTerm st : sts)
        {
            EXPECT_FALSE(st.getProcessVarName().empty());
            EXPECT_FALSE(st.getParamName().empty());
            EXPECT_FALSE(DataHolderLib::ConditionType::NONE == st.getType());
            EXPECT_FALSE(st.getBaseObjName().empty());
            if (st.getBaseObjType() == DataHolderLib::BaseObjType::GEOMETRY)
                EXPECT_FALSE(st.getObjName().empty());
        }
    }
}
