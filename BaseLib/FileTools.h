/**
 * \file
 * \author Lars Bilke
 * \date   Apr. 2010
 * \brief Filename manipulation routines.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "Logging.h"

namespace BaseLib
{
/**
 * \brief Returns true if given file exists.
 *
 * \param strFilename         the file name
 */
bool IsFileExisting(const std::string &strFilename);

/**
 * Returns the begin and end position of the string enclosed in open_char and
 * close_char and the enclosed string itself. Search starts at position pos
 * within the string str. Nested open_char and close_char are not handled
 * correctly.
 */
std::tuple<std::string, std::string::size_type, std::string::size_type>
getParenthesizedString(std::string const& in,
                       char const open_char,
                       char const close_char,
                       std::string::size_type pos);

std::string constructFormattedFileName(std::string const& format_specification,
                                       std::string const& mesh_name,
                                       int const timestep,
                                       double const t,
                                       int const iteration);

/**
 * \brief write value as binary into the given output stream
 *
 * \tparam T    data type of the value
 * \param out   output stream, have to be opened in binary mode
 * \param val   value
 */
template <typename T> void writeValueBinary(std::ostream &out, T const& val)
{
    out.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template <typename T>
T swapEndianness(T const& v)
{
    union
    {
        T v;
        char c[sizeof(T)];
    } a, b;

    a.v = v;
    for (unsigned short i = 0; i < sizeof(T); i++)
    {
        b.c[i] = a.c[sizeof(T) - i - 1];
    }

    return b.v;
}

double swapEndianness(double const& v);


template <typename T>
T readBinaryValue(std::istream& in)
{
    T v;
    in.read(reinterpret_cast<char*>(&v), sizeof(T));
    return v;
}

template <typename T>
std::vector<T> readBinaryArray(std::string const& filename, std::size_t const n)
{
    std::ifstream in(filename.c_str());
    if (!in) {
        ERR("readBinaryArray(): Error while reading from file '{:s}'.",
            filename);
        ERR("Could not open file '{:s}' for input.", filename);
        in.close();
        return std::vector<T>();
    }

    std::vector<T> result;
    result.reserve(n);

    for (std::size_t p = 0; in && !in.eof() && p < n; ++p)
    {
        result.push_back(BaseLib::readBinaryValue<T>(in));
    }

    if (result.size() == n)
    {
        return result;
    }

    ERR("readBinaryArray(): Error while reading from file '{:s}'.", filename);
    ERR("Read different number of values. Expected {:d}, got {:d}.",
        n,
        result.size());

    if (!in.eof())
    {
        ERR("EOF reached.\n");
    }

    return std::vector<T>();
}

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
std::string getFileExtension(std::string const& path);

/**
 * Compares filename's extension with query extension. The comparison is case
 * insensitive.
 */
bool hasFileExtension(std::string const& extension,
                      std::string const& filename);

/**
 * Checks if file_name already contains a qualified path and if not copies the
 * path from source.
 */
std::string copyPathToFileName(const std::string &file_name,
                               const std::string &source);

/** Returns a string with file extension as found by getFileExtension()
 * dropped.
 */
std::string dropFileExtension(std::string const& filename);

/**
 * Extracts the path of a pathname.
 *
 * Returns a string up to the last path separator not including it.
 */
std::string extractPath(std::string const& pathname);

/**
 * Concat two paths. Does not check for validity.
 */
std::string joinPaths(std::string const& pathA, std::string const& pathB);

/// Returns the directory where the prj file resides.
std::string const& getProjectDirectory();

/// Sets the project directory.
void setProjectDirectory(std::string const& dir);

/// Removes a file. If a file does not exist nothing will happen, other errors
/// lead to OGS_FATAL call.
void removeFile(std::string const& filename);

/// Remove files. If a file does not exist nothing will happen, other errors
/// lead to OGS_FATAL call.
void removeFiles(std::vector<std::string> const& files);
}  // end namespace BaseLib
