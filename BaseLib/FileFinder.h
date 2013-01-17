/**
 * \file
 * \author Karsten Rink
 * \date   2010-10-26
 * \brief  Definition of the FileFinder class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef FILEFINDER_H
#define FILEFINDER_H

#include <fstream>
#include <string>
#include <vector>

// ThirdParty/logog
#include "logog/include/logog.hpp"

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
	/// Constructor
	FileFinder() {}

	/**
	 * \brief Adds another directory to the search-space.
	 * If the given directory does not end with a slash one will be appended.
	 */
	void addDirectory(std::string const& dir)
	{
		if (dir[dir.size() - 1] != '/')
			_directories.push_back(std::string(dir + "/"));
		else
			_directories.push_back(dir);
	}

	/**
	 * Given a filename, this method will return the complete path where this file can be found.
	 * If the file is located in more than one of the directories in the search list, only the
	 * first location will be returned.
	 */
	std::string getPath(std::string const& filename) const
	{
		if (_directories.empty())
			ERR("FileFinder::getPath(): No directories set.");

		for (std::vector<std::string>::const_iterator it = _directories.begin(); it
		     != _directories.end(); ++it)
		{
			std::string testDir(*it);
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

private:
	std::vector<std::string> _directories;
};
} // end namespace BaseLib

#endif // FILEFINDER_H
