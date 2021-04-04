/**
 * \file
 * \author
 * \date
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include "BaseLib/FileTools.h"

#ifdef WIN32
TEST(BaseLib, CopyPathToFileNameWin)
{
    ASSERT_EQ("extend\\file", BaseLib::copyPathToFileName("file", "extend"));
    ASSERT_EQ("path\\file",
              BaseLib::copyPathToFileName("path\\file", "extend"));
    ASSERT_EQ("extend\\file", BaseLib::copyPathToFileName("file", "extend\\"));
    ASSERT_EQ("path\\file",
              BaseLib::copyPathToFileName("path\\file", "extend\\"));
    ASSERT_EQ("extend\\smth\\file",
              BaseLib::copyPathToFileName("file", "extend\\smth"));
    ASSERT_EQ("path\\file",
              BaseLib::copyPathToFileName("path\\file", "extend\\smth"));
}
#else
TEST(BaseLib, CopyPathToFileNameUnix)
{
    ASSERT_EQ("extend/file", BaseLib::copyPathToFileName("file", "extend"));
    ASSERT_EQ("path/file", BaseLib::copyPathToFileName("path/file", "extend"));
    ASSERT_EQ("extend/file", BaseLib::copyPathToFileName("file", "extend/"));
    ASSERT_EQ("path/file", BaseLib::copyPathToFileName("path/file", "extend/"));

    ASSERT_EQ("extend/smth/file",
              BaseLib::copyPathToFileName("file", "extend/smth"));
    ASSERT_EQ("path/file",
              BaseLib::copyPathToFileName("path/file", "extend/smth"));
}
#endif
