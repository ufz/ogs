/**
 * \file
 * \author Lars Bilke
 * \date   Apr. 2010
 * \brief Filename manipulation routines.
 *
 * \copyright
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef FILETOOLS_H
#define FILETOOLS_H

#include <fstream>
#include <string>

namespace BaseLib
{
/**
 * \brief Returns true if given file exists. From http://www.techbytes.ca/techbyte103.html
 *
 * \param strFilename         the file name
 */
bool IsFileExisting(const std::string &strFilename);

/**
 * \brief write value as binary into the given output stream
 *
 * \param T    data type of the value
 * \param out   output stream, have to be opened in binary mode
 * \param val   value
 */
template <typename T> void writeValueBinary(std::ostream &out, T const& val)
{
	out.write((const char*)&val, sizeof(T));
}

/**
 * \brief truncate a file
 *
 * \param file_path         the file name
 */
void truncateFile( std::string const& file_path);

/**
 * Extracts basename from given pathname with extension.
 *
 * Returns a string containing everything after the last path separator.
 * If the the pathname does not contain a path separator original pathname is
 * returned.
 */
std::string extractBaseName(std::string const& pathname);

/**
 * Extracts basename from given pathname without its extension.
 *
 *  Same as extractBaseName(), but drops the file extension too.
 */
std::string extractBaseNameWithoutExtension(std::string const& pathname);

/**
 * Extract extension from filename
 */
std::string getFileExtension(std::string const& filename);

/**
 * Compares filename's extension with query extension. The comparison is case
 * insensitive done by converting to upper case with the std::toupper()
 * function.
 */
bool hasFileExtension(std::string const& extension, std::string const& filename);

/** Returns a string with file extension as found by getFileExtension()
 * dropped.
 */
std::string dropFileExtension(std::string const& filename);

/**
 * Checks if file_name already contains a qualified path and if not copies the
 * path from source.
 */
std::string copyPathToFileName(const std::string &file_name, const std::string &source);

/**
 * Extracts the path of a pathname.
 *
 * Returns a string up to the last path separator not including it.
 */
std::string extractPath(std::string const& pathname);
} // end namespace BaseLib

#endif // FILETOOLS_H
