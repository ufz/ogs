/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "BaseLib/IO/readStringListFromFile.h"
#include "InfoLib/TestInfo.h"

TEST(BaseLibReadStringListFromFile, testreader)
{
    std::string const filename = TestInfoLib::TestInfo::data_path +
                                 "/MeshGeoToolsLib/Ammer/AmmerLayers.txt";
    auto const list = BaseLib::IO::readStringListFromFile(filename);
    ASSERT_EQ(4, list.size());
    for (std::size_t i = 0; i < list.size(); ++i)
    {
        ASSERT_EQ(std::string("layer" + std::to_string(i) + ".vtu"), list[i]);
    }
}
