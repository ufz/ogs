/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FileFinder.h
 *
 * Created on 2010-10-26 by Karsten Rink
 */

#ifndef FILEFINDER_H
#define FILEFINDER_H

#include <iostream>
#include <fstream>
#include <string>
#include <list>

namespace BaseLib {
/**
 * FileFinder stores a list of directories and will return the complete path
 * for a given filename if the corresponding file is found in any of these
 * directories.
 */
class FileFinder
{
public:
	/// Constructor
	FileFinder() {};

	/**
	 * \brief Adds another directory to the search-space.
	 * If the given directory does not end with a slash one will be appended.
	 */
	void addDirectory(std::string dir)
	{
		if (dir[dir.size()-1] != '/') dir.append("/");
		_directories.push_back(dir);
	};

	/**
	 * Given a filename, this method will return the complete path where this file can be found.
	 * If the file is located in more than one of the directories in the search list, only the
	 * first location will be returned.
	 */
	std::string getPath(std::string filename) const
	{
		if (_directories.empty()) std::cout << "Error: FileFinder::getPath() -- directory list is empty." << std::endl;
		for (std::list<std::string>::const_iterator it = _directories.begin(); it != _directories.end(); ++it)
		{
			std::string testDir(*it);
			std::ifstream is(testDir.append(filename).c_str());
			if (is.good()) return testDir;
		}
		std::cout << "Error: FileFinder::getPath() -- file not found." << std::endl;
		return filename;
	};

private:
	std::list<std::string> _directories;
};

} // end namespace BaseLib

#endif // FILEFINDER_H
