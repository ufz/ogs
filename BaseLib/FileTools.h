/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file FileTools.h
 *
 * Created on 2010-04-26 by Lars Bilke
 *
 */


#ifndef FILETOOLS_H
#define FILETOOLS_H

#include <string>
#include <fstream>

namespace BaseLib
{

/** Finds the position of last file path separator.
 * Checks for unix or windows file path separators in given path and returns the
 * position of the last one or std::string::npos if no file path separator was
 * found.
 */
size_t findLastPathSeparator(std::string const& path);

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
 * Extract the filename from a path
 */
std::string getFileNameFromPath(const std::string &str, bool with_extension = false);

/**
 * Extract the file type / suffix from a path
 */
std::string getSuffixFromPath(const std::string &str);


/**
 * Checks if file_name already contains a qualified path and if not copies the path from source.
 */
std::string copyPathToFileName(const std::string &file_name, const std::string &source);

/**
 * Extracts the path of a pathname.
 */
std::string extractPath(std::string const& pathname);
} // end namespace BaseLib

#endif // FILETOOLS_H
