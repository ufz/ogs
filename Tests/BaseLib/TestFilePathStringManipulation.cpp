/**
 * \file
 * \author
 * \date
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "gtest/gtest.h"

#include "BaseLib/FileTools.h"

#ifdef WIN32
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
#else
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
#endif

#ifdef WIN32
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
#else
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
#endif
