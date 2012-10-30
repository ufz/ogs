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

/**
 * \brief Returns true if given file exists. From http://www.techbytes.ca/techbyte103.html
 *
 * \param strFilename         the file name
 */
bool IsFileExisting(const std::string &strFilename);

/**
 * \brief write value as binary into the given output stream
 *
 * \tparam T    data type of the value
 * \param out   output stream
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

} // end namespace BaseLib

#endif // FILETOOLS_H
