/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include "BaseLib/FileTools.h"

#ifdef WIN32
TEST(BaseLib, JoinPathsWin)
{
    ASSERT_EQ("extend\\file", BaseLib::joinPaths("extend", "file"));
    ASSERT_EQ("extend\\path\\file", BaseLib::joinPaths("extend", "path\\file"));
    ASSERT_EQ("extend\\file", BaseLib::joinPaths("extend\\", "file"));
    ASSERT_EQ("extend\\path\\file",
              BaseLib::joinPaths("extend\\", "path\\file"));
    ASSERT_EQ("extend\\smth\\file", BaseLib::joinPaths("extend\\smth", "file"));
    ASSERT_EQ("extend\\smth\\path\\file",
              BaseLib::joinPaths("extend\\smth", "path\\file"));
}
#else
TEST(BaseLib, JoinPathsUnix)
{
    ASSERT_EQ("extend/file", BaseLib::joinPaths("extend", "file"));
    ASSERT_EQ("extend/path/file", BaseLib::joinPaths("extend", "path/file"));
    ASSERT_EQ("extend/file", BaseLib::joinPaths("extend/", "file"));
    ASSERT_EQ("extend/path/file", BaseLib::joinPaths("extend/", "path/file"));

    ASSERT_EQ("extend/smth/file", BaseLib::joinPaths("extend/smth", "file"));
    ASSERT_EQ("extend/smth/path/file",
              BaseLib::joinPaths("extend/smth", "path/file"));
}
#endif
