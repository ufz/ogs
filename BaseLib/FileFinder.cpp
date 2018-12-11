/**
 * \file
 * \author Karsten Rink
 * \date   2010-10-26
 * \brief  Definition of the FileFinder class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FileFinder.h"

#include <fstream>

#include "Logging.h"

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
    {
        addDirectory(dir);
    }
}

void FileFinder::addDirectory(std::string const& dir)
{
    if (dir.empty())
    {
        return;
    }

    if (dir[dir.size() - 1] != '/')
    {
        _directories.emplace_back(dir + "/");
    }
    else
    {
        _directories.push_back(dir);
    }
}
}  // namespace BaseLib
