// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
