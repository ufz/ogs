/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "gtest/gtest.h"

#include "FileTools.h"
#include "FileTools.cpp"

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

TEST(BaseLib, DropFileExtensionWin)
{
	ASSERT_EQ ( BaseLib::dropFileExtension("file"), "file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("\\file"), "\\file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("path\\"), "path\\" );
	ASSERT_EQ ( BaseLib::dropFileExtension("\\path\\"), "\\path\\" );
	ASSERT_EQ ( BaseLib::dropFileExtension("path\\file"), "path\\file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("\\path\\file"), "\\path\\file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("\\path\\path\\file"), "\\path\\path\\file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("\\path\\path\\path\\"), "\\path\\path\\path\\" );

	ASSERT_EQ ( BaseLib::dropFileExtension("file.ext"), "file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("\\file.ext"), "\\file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("path.ext\\"), "path.ext\\" );
	ASSERT_EQ ( BaseLib::dropFileExtension("\\path.ext\\"), "\\path.ext\\" );
	ASSERT_EQ ( BaseLib::dropFileExtension("path\\file.ext"), "path\\file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("\\path\\file.ext"), "\\path\\file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("\\path\\path\\file.ext"), "\\path\\path\\file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("\\path\\path\\path.ext\\"), "\\path\\path\\path.ext\\" );

	ASSERT_EQ ( BaseLib::dropFileExtension("path.wrong\\file.ext"), "path.wrong\\file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("\\path.wrong\\file.ext"), "\\path.wrong\\file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("\\path.wrong0\\path.wrong\\file.ext"), "\\path.wrong0\\path.wrong\\file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("\\path.wrong0\\path.wrong\\path.ext\\"), "\\path.wrong0\\path.wrong\\path.ext\\" );
}

TEST(BaseLib, DropFileExtensionUnix)
{
	ASSERT_EQ ( BaseLib::dropFileExtension("file"), "file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("/file"), "/file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("path/"), "path/" );
	ASSERT_EQ ( BaseLib::dropFileExtension("/path/"), "/path/" );
	ASSERT_EQ ( BaseLib::dropFileExtension("path/file"), "path/file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("/path/file"), "/path/file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("/path/path/file"), "/path/path/file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("/path/path/path/"), "/path/path/path/" );

	ASSERT_EQ ( BaseLib::dropFileExtension("file.ext"), "file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("/file.ext"), "/file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("path.ext/"), "path.ext/" );
	ASSERT_EQ ( BaseLib::dropFileExtension("/path.ext/"), "/path.ext/" );
	ASSERT_EQ ( BaseLib::dropFileExtension("path/file.ext"), "path/file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("/path/file.ext"), "/path/file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("/path/path/file.ext"), "/path/path/file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("/path/path/path.ext/"), "/path/path/path.ext/" );

	ASSERT_EQ ( BaseLib::dropFileExtension("path.wrong/file.ext"), "path.wrong/file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("/path.wrong/file.ext"), "/path.wrong/file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("/path.wrong0/path.wrong/file.ext"), "/path.wrong0/path.wrong/file" );
	ASSERT_EQ ( BaseLib::dropFileExtension("/path.wrong0/path.wrong/path.ext/"), "/path.wrong0/path.wrong/path.ext/" );
}

TEST(BaseLib, GetFileExtensionWin)
{
	ASSERT_EQ ( BaseLib::getFileExtension("file"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("\\file"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("path\\"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("\\path\\"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("path\\file"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("\\path\\file"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("\\path\\path\\file"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("\\path\\path\\path\\"), "" );

	ASSERT_EQ ( BaseLib::getFileExtension("file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("path.ext\\"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("\\path.ext\\"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("path\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("\\path\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("\\path\\path\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("\\path\\path\\path.ext\\"), "" );

	ASSERT_EQ ( BaseLib::getFileExtension("path.wrong\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("\\path.wrong\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("\\path.wrong0\\path.wrong\\file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("\\path.wrong0\\path.wrong\\path.ext\\"), "" );
}

TEST(BaseLib, getFileExtensionUnix)
{
	ASSERT_EQ ( BaseLib::getFileExtension("file"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("/file"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("path/"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("/path/"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("path/file"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("/path/file"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("/path/path/file"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("/path/path/path/"), "" );

	ASSERT_EQ ( BaseLib::getFileExtension("file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("path.ext/"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("/path.ext/"), "" );
	ASSERT_EQ ( BaseLib::getFileExtension("path/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("/path/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("/path/path/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("/path/path/path.ext/"), "" );

	ASSERT_EQ ( BaseLib::getFileExtension("path.wrong/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("/path.wrong/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("/path.wrong0/path.wrong/file.ext"), "ext" );
	ASSERT_EQ ( BaseLib::getFileExtension("/path.wrong0/path.wrong/path.ext/"), "" );
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

TEST(BaseLib, ExtractBaseNameWithoutExtensionWin)
{
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("\\file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("path\\"), "" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("\\path\\"), "" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("path\\file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("\\path\\file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("\\path\\path\\file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("\\path\\path\\path\\"), "" );

	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("\\file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("path.ext\\"), "" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("\\path.ext\\"), "" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("path\\file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("\\path\\file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("\\path\\path\\file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("\\path\\path\\path.ext\\"), "" );

	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("path.wrong\\file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("\\path.wrong\\file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("\\path.wrong0\\path.wrong\\file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("\\path.wrong0\\path.wrong\\path.ext\\"), "" );
}

TEST(BaseLib, ExtractBaseNameWithoutExtensionUnix)
{
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("/file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("path/"), "" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("/path/"), "" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("path/file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("/path/file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("/path/path/file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("/path/path/path/"), "" );

	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("/file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("path.ext/"), "" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("/path.ext/"), "" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("path/file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("/path/file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("/path/path/file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("/path/path/path.ext/"), "" );

	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("path.wrong/file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("/path.wrong/file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("/path.wrong0/path.wrong/file.ext"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseNameWithoutExtension("/path.wrong0/path.wrong/path.ext/"), "" );
}

TEST(BaseLib, ExtractBaseNameWin)
{
	ASSERT_EQ ( BaseLib::extractBaseName("file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseName("\\file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseName("path\\"), "" );
	ASSERT_EQ ( BaseLib::extractBaseName("\\path\\"), "" );
	ASSERT_EQ ( BaseLib::extractBaseName("path\\file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseName("\\path\\file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseName("\\path\\path\\file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseName("\\path\\path\\path\\"), "" );

	ASSERT_EQ ( BaseLib::extractBaseName("file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("\\file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("path.ext\\"), "" );
	ASSERT_EQ ( BaseLib::extractBaseName("\\path.ext\\"), "" );
	ASSERT_EQ ( BaseLib::extractBaseName("path\\file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("\\path\\file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("\\path\\path\\file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("\\path\\path\\path.ext\\"), "" );

	ASSERT_EQ ( BaseLib::extractBaseName("path.wrong\\file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("\\path.wrong\\file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("\\path.wrong0\\path.wrong\\file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("\\path.wrong0\\path.wrong\\path.ext\\"), "" );
}

TEST(BaseLib, HasFileExtensionWin)
{
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "file"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "\\file"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "path\\"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "\\path\\"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "path\\file"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "\\path\\file"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "\\path\\path\\file"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "\\path\\path\\path\\"));

	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "\\file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "path.ext\\"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "\\path.ext\\"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "path\\file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "\\path\\file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "\\path\\path\\file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "\\path\\path\\path.ext\\"));

	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "path.wrong\\file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "\\path.wrong\\file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "\\path.wrong0\\path.wrong\\file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "\\path.wrong0\\path.wrong\\path.ext\\"));

	ASSERT_TRUE ( BaseLib::hasFileExtension("EXT", "file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("EXT", "file.EXT"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "file.EXT"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("Ext", "file.exT"));

	ASSERT_TRUE ( BaseLib::hasFileExtension("EXT", "path\\file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("EXT", "path\\file.EXT"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "path\\file.EXT"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("Ext", "path\\file.exT"));
}

TEST(BaseLib, HasFileExtensionUnix)
{
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "file"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "/file"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "path/"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "/path/"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "path/file"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "/path/file"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "/path/path/file"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "/path/path/path/"));

	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "/file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "path.ext/"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "/path.ext/"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "path/file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "/path/file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "/path/path/file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "/path/path/path.ext/"));

	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "path.wrong/file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "/path.wrong/file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "/path.wrong0/path.wrong/file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("", "/path.wrong0/path.wrong/path.ext/"));

	ASSERT_TRUE ( BaseLib::hasFileExtension("EXT", "file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("EXT", "file.EXT"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "file.EXT"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("Ext", "file.exT"));

	ASSERT_TRUE ( BaseLib::hasFileExtension("EXT", "path/file.ext"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("EXT", "path/file.EXT"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("ext", "path/file.EXT"));
	ASSERT_TRUE ( BaseLib::hasFileExtension("Ext", "path/file.exT"));
}

TEST(BaseLib, ExtractBaseNameUnix)
{
	ASSERT_EQ ( BaseLib::extractBaseName("file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseName("/file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseName("path/"), "" );
	ASSERT_EQ ( BaseLib::extractBaseName("/path/"), "" );
	ASSERT_EQ ( BaseLib::extractBaseName("path/file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseName("/path/file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseName("/path/path/file"), "file" );
	ASSERT_EQ ( BaseLib::extractBaseName("/path/path/path/"), "" );

	ASSERT_EQ ( BaseLib::extractBaseName("file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("/file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("path.ext/"), "" );
	ASSERT_EQ ( BaseLib::extractBaseName("/path.ext/"), "" );
	ASSERT_EQ ( BaseLib::extractBaseName("path/file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("/path/file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("/path/path/file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("/path/path/path.ext/"), "" );

	ASSERT_EQ ( BaseLib::extractBaseName("path.wrong/file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("/path.wrong/file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("/path.wrong0/path.wrong/file.ext"), "file.ext" );
	ASSERT_EQ ( BaseLib::extractBaseName("/path.wrong0/path.wrong/path.ext/"), "" );
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
