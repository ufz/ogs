// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <iosfwd>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

namespace BaseLib
{
/**
 * \brief Returns true if given file exists.
 *
 * \param strFilename         the file name
 */
bool IsFileExisting(const std::string& strFilename);

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
                                       int const timestep, double const t,
                                       int const iteration,
                                       bool const converged);

std::vector<double> readDoublesFromBinaryFile(
    const std::string& filename, std::string const& project_directory);

/**
 * \brief write value as binary into the given output stream
 *
 * \tparam T    data type of the value
 * \param out   output stream, have to be opened in binary mode
 * \param val   value
 */
template <typename T>
void writeValueBinary(std::ostream& out, T const& val);

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
T readBinaryValue(std::istream& in);

extern template float readBinaryValue<float>(std::istream&);
extern template double readBinaryValue<double>(std::istream&);

template <typename T>
std::vector<T> readBinaryVector(
    std::string const& filename, std::size_t const start_element = 0,
    std::size_t const num_elements = std::numeric_limits<std::size_t>::max());

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

/// Removes a file. If a file does not exist nothing will happen, other errors
/// lead to OGS_FATAL call.
void removeFile(std::string const& filename);

/// Remove files. If a file does not exist nothing will happen, other errors
/// lead to OGS_FATAL call.
void removeFiles(std::vector<std::string> const& files);

// Creates the given directory. Returns true on success.
bool createOutputDirectory(std::string const& dir);
}  // end namespace BaseLib
