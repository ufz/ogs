// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
