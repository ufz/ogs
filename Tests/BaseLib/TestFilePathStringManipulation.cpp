/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "gtest.h"

#include "FileTools.h"

TEST(BaseLib, FindLastPathSeparatorWin)
{
	ASSERT_EQ ( BaseLib::findLastPathSeparator("file"), std::string::npos );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("\\file"), 0 );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("path\\"), 4 );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("\\path\\"), 5 );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("path\\file"), 4 );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("\\path\\file"), 5 );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("\\path\\path\\file"), 10 );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("\\path\\path\\path\\"), 15 );
}

TEST(BaseLib, FindLastPathSeparatorUnix)
{
	ASSERT_EQ ( BaseLib::findLastPathSeparator("file"), std::string::npos );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("/file"), 0 );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("path/"), 4 );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("/path/"), 5 );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("path/file"), 4 );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("/path/file"), 5 );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("/path/path/file"), 10 );
	ASSERT_EQ ( BaseLib::findLastPathSeparator("/path/path/path/"), 15 );
}

TEST(BaseLib, GetFileNameFromPathWin)
{
	ASSERT_EQ ( BaseLib::getFileNameFromPath("file", true), "file" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("\\file", true), "file" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("path\\", true), "" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("\\path\\", true), "" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("path\\file", true), "file" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("\\path\\file", true), "file" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("\\path\\path\\file", true), "file" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("\\path\\path\\path\\", true), "" );

	ASSERT_EQ ( BaseLib::getFileNameFromPath("file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("\\file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("path.ext\\", true), "" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("\\path.ext\\", true), "" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("path\\file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("\\path\\file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("\\path\\path\\file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("\\path\\path\\path.ext\\", true), "" );

	ASSERT_EQ ( BaseLib::getFileNameFromPath("path.wrong\\file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("\\path.wrong\\file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("\\path.wrong0\\path.wrong\\file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("\\path.wrong0\\path.wrong\\path.ext\\", true), "" );
}

TEST(BaseLib, GetFileNameFromPathUnix)
{
	ASSERT_EQ ( BaseLib::getFileNameFromPath("file", true), "file" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("/file", true), "file" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("path/", true), "" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("/path/", true), "" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("path/file", true), "file" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("/path/file", true), "file" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("/path/path/file", true), "file" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("/path/path/path/", true), "" );

	ASSERT_EQ ( BaseLib::getFileNameFromPath("file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("/file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("path.ext/", true), "" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("/path.ext/", true), "" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("path/file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("/path/file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("/path/path/file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("/path/path/path.ext/", true), "" );

	ASSERT_EQ ( BaseLib::getFileNameFromPath("path.wrong/file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("/path.wrong/file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("/path.wrong0/path.wrong/file.ext", true), "file.ext" );
	ASSERT_EQ ( BaseLib::getFileNameFromPath("/path.wrong0/path.wrong/path.ext/", true), "" );
}

TEST(BaseLib, GetSuffixFromPathWin)
{
	ASSERT_EQ ( BaseLib::getSuffixFromPath("file"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("\\file"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("path\\"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("\\path\\"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("path\\file"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("\\path\\file"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("\\path\\path\\file"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("\\path\\path\\path\\"), "" );

	ASSERT_EQ ( BaseLib::getSuffixFromPath("file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("path.ext\\"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("\\path.ext\\"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("path\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("\\path\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("\\path\\path\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("\\path\\path\\path.ext\\"), "" );

	ASSERT_EQ ( BaseLib::getSuffixFromPath("path.wrong\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("\\path.wrong\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("\\path.wrong0\\path.wrong\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("\\path.wrong0\\path.wrong\\path.ext\\"), "" );
}

TEST(BaseLib, GetSuffixFromPathUnix)
{
	ASSERT_EQ ( BaseLib::getSuffixFromPath("file"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("/file"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("path/"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("/path/"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("path/file"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("/path/file"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("/path/path/file"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("/path/path/path/"), "" );

	ASSERT_EQ ( BaseLib::getSuffixFromPath("file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("path.ext/"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("/path.ext/"), "" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("path/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("/path/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("/path/path/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("/path/path/path.ext/"), "" );

	ASSERT_EQ ( BaseLib::getSuffixFromPath("path.wrong/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("/path.wrong/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("/path.wrong0/path.wrong/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getSuffixFromPath("/path.wrong0/path.wrong/path.ext/"), "" );
}

TEST(BaseLib, CopyPathToFileNameWin)
{
	ASSERT_EQ ( BaseLib::copyPathToFileName("file", "extend"), "file" );
	ASSERT_EQ ( BaseLib::copyPathToFileName("path\\file", "extend"), "path\\file" );

	ASSERT_EQ ( BaseLib::copyPathToFileName("file", "extend\\"), "extend\\file" );
	ASSERT_EQ ( BaseLib::copyPathToFileName("path\\file", "extend\\"), "path\\file" );

	ASSERT_EQ ( BaseLib::copyPathToFileName("file", "extend\\smth"), "extend\\file" );
	ASSERT_EQ ( BaseLib::copyPathToFileName("path\\file", "extend\\smth"), "path\\file" );
}

TEST(BaseLib, CopyPathToFileNameUnix)
{
	ASSERT_EQ ( BaseLib::copyPathToFileName("file", "extend"), "file" );
	ASSERT_EQ ( BaseLib::copyPathToFileName("path/file", "extend"), "path/file" );

	ASSERT_EQ ( BaseLib::copyPathToFileName("file", "extend/"), "extend/file" );
	ASSERT_EQ ( BaseLib::copyPathToFileName("path/file", "extend/"), "path/file" );

	ASSERT_EQ ( BaseLib::copyPathToFileName("file", "extend/smth"), "extend/file" );
	ASSERT_EQ ( BaseLib::copyPathToFileName("path/file", "extend/smth"), "path/file" );
}

TEST(BaseLib, ExtractPathWin)
{
	ASSERT_EQ ( BaseLib::extractPath("file"), "" );
	ASSERT_EQ ( BaseLib::extractPath("/file"), "/" );
	ASSERT_EQ ( BaseLib::extractPath("path/"), "path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/"), "/path/" );
	ASSERT_EQ ( BaseLib::extractPath("path/file"), "path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/file"), "/path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/path/file"), "/path/path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/path/path/"), "/path/path/path/" );

	ASSERT_EQ ( BaseLib::extractPath("file.ext"), "" );
	ASSERT_EQ ( BaseLib::extractPath("/file.ext"), "/" );
	ASSERT_EQ ( BaseLib::extractPath("path.ext/"), "path.ext/" );
	ASSERT_EQ ( BaseLib::extractPath("/path.ext/"), "/path.ext/" );
	ASSERT_EQ ( BaseLib::extractPath("path/file.ext"), "path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/file.ext"), "/path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/path/file.ext"), "/path/path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/path/path.ext/"), "/path/path/path.ext/" );
}

TEST(BaseLib, ExtractPathUnix)
{
	ASSERT_EQ ( BaseLib::extractPath("file"), "" );
	ASSERT_EQ ( BaseLib::extractPath("/file"), "/" );
	ASSERT_EQ ( BaseLib::extractPath("path/"), "path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/"), "/path/" );
	ASSERT_EQ ( BaseLib::extractPath("path/file"), "path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/file"), "/path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/path/file"), "/path/path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/path/path/"), "/path/path/path/" );

	ASSERT_EQ ( BaseLib::extractPath("file.ext"), "" );
	ASSERT_EQ ( BaseLib::extractPath("/file.ext"), "/" );
	ASSERT_EQ ( BaseLib::extractPath("path.ext/"), "path.ext/" );
	ASSERT_EQ ( BaseLib::extractPath("/path.ext/"), "/path.ext/" );
	ASSERT_EQ ( BaseLib::extractPath("path/file.ext"), "path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/file.ext"), "/path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/path/file.ext"), "/path/path/" );
	ASSERT_EQ ( BaseLib::extractPath("/path/path/path.ext/"), "/path/path/path.ext/" );
}
