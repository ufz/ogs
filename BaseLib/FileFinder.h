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

#pragma once

#include <initializer_list>
#include <string>
#include <vector>


namespace BaseLib
{
/**
 * FileFinder stores a list of directories and will return the complete path
 * for a given filename if the corresponding file is found in any of these
 * directories.
 */
class FileFinder
{
public:
    /// Constructor having current directory (.) as the search-space
    FileFinder();

    /**
     * Construct with the given directory paths in addition to current directory (.)
     *
     * @param dirs   an initializer list of additional directory paths to the search-space
     */
    FileFinder(std::initializer_list<std::string> dirs);

    /**
     * \brief Adds another directory to the search-space.
     * If the given directory does not end with a slash one will be appended.
     */
    void addDirectory(std::string const& dir);

    /**
     * Given a filename, this method will return the complete path where this file can be found.
     * If the file is located in more than one of the directories in the search list, only the
     * first location will be returned.
     */
    std::string getPath(std::string const& filename) const;

private:
    std::vector<std::string> _directories;
};
} // end namespace BaseLib
