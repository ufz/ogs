/**
 * \file
 * \author Karsten Rink
 * \date   2010-10-26
 * \brief  Definition of the FileFinder class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FileFinder.h"

#include <fstream>

#include <logog/include/logog.hpp>


namespace BaseLib
{

FileFinder::FileFinder()
{
    addDirectory(".");
}

FileFinder::FileFinder(std::initializer_list<std::string> dirs)
{
    addDirectory(".");
    for (auto const& dir : dirs)
        addDirectory(dir);
}

void FileFinder::addDirectory(std::string const& dir)
{
    if (dir.empty())
        return;

    if (dir[dir.size() - 1] != '/')
        _directories.emplace_back(dir + "/");
    else
        _directories.push_back(dir);
}

std::string FileFinder::getPath(std::string const& filename) const
{
    if (_directories.empty())
        ERR("FileFinder::getPath(): No directories set.");

    for (auto const& dir : _directories)
    {
        std::string testDir(dir);
        std::ifstream is(testDir.append(filename).c_str());
        if (is.good())
        {
            is.close();
            return testDir;
        }
    }
    ERR("FileFinder::getPath(): File not found.");
    return filename;
}

} // end namespace BaseLib
