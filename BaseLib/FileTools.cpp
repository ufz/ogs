/**
 * \file
 * \author Lars Bilke
 * \date   Apr. 2010
 * \brief Filename manipulation routines implementation.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FileTools.h"
#include "StringTools.h"

#include <sys/stat.h>
#include <boost/filesystem.hpp>

namespace BaseLib
{
/**
 * Returns true if given file exists. From http://www.techbytes.ca/techbyte103.html
 */
bool IsFileExisting(const std::string &strFilename)
{
	if(boost::filesystem::exists(strFilename))
		return true;
	else
		return false;
}

/**
 * \brief truncate a file
 */
void truncateFile( std::string const& filename)
{
	std::ofstream ofs(filename.c_str(), std::ios_base::trunc);
	ofs.close();
}

/** Finds the position of last file path separator.
 * Checks for unix or windows file path separators in given path and returns the
 * position of the last one or std::string::npos if no file path separator was
 * found.
 */
static
size_t findLastPathSeparator(std::string const& path)
{
	return path.find_last_of("/\\");
}

/** Finds the position of last dot.
 * This could be used to extract file extension.
 */
static
size_t findLastDot(std::string const& path)
{
	return path.find_last_of(".");
}

std::string dropFileExtension(std::string const& filename)
{
	// Look for dots in filename.
	const size_t p = findLastDot(filename);
	if (p == std::string::npos)
		return filename;

	// Check position of the last path separator.
	const size_t s = findLastPathSeparator(filename);
	if (s != std::string::npos && p < s)
		return filename;

	return filename.substr(0, p);
}

std::string extractBaseName(std::string const& pathname)
{
	const size_t p = findLastPathSeparator(pathname);
	return pathname.substr(p + 1);
}

std::string extractBaseNameWithoutExtension(std::string const& pathname)
{
	std::string basename = extractBaseName(pathname);
	return dropFileExtension(basename);
}

std::string getFileExtension(const std::string &path)
{
	const std::string str = extractBaseName(path);
	const size_t p = findLastDot(str);
	if (p == std::string::npos)
		return std::string();
	return str.substr(p + 1);
}

bool hasFileExtension(std::string const& extension, std::string const& filename)
{
	std::string ext = stringToUpper(extension); // Copy for modification.
	std::string file_ext = stringToUpper(getFileExtension(filename));

	return ext == file_ext;
}

std::string copyPathToFileName(const std::string &file_name,
                               const std::string &source)
{
	// check if file_name already contains a full path
	const size_t pos = findLastPathSeparator(file_name);
	if (pos != std::string::npos)
		return file_name;

	return BaseLib::extractPath(source).append(file_name);
}

std::string extractPath(std::string const& pathname)
{
	const size_t pos = findLastPathSeparator(pathname);
	return pathname.substr(0, pos + 1);
}
} // end namespace BaseLib
