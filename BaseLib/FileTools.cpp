/**
 * \file
 * \author Lars Bilke
 * \date   Apr. 2010
 * \brief Filename manipulation routines implementation.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FileTools.h"
#include "StringTools.h"

#include <sys/stat.h>
#include <boost/algorithm/string.hpp>

namespace BaseLib
{
/**
 * Returns true if given file exists. From http://www.techbytes.ca/techbyte103.html
 */
bool IsFileExisting(const std::string &strFilename)
{
    struct stat buffer;
    return (stat (strFilename.c_str(), &buffer) == 0);
}

double swapEndianness(double const& v)
{
    union
    {
        double v;
        char c[sizeof(double)];
    } a, b;

    a.v = v;
    for (unsigned short i = 0; i < sizeof(double)/2; i++)
        b.c[i] = a.c[sizeof(double)/2 - i - 1];

    for (unsigned short i = sizeof(double)/2; i < sizeof(double); i++)
        b.c[i] = a.c[sizeof(double)+sizeof(double)/2 - i - 1];

    return b.v;
}

/**
 * \brief truncate a file
 */
void truncateFile( std::string const& filename)
{
    std::ofstream ofs(filename.c_str(), std::ios_base::trunc);
    ofs.close();
}

namespace
{

/** Finds the position of last file path separator.
 * Checks for unix or windows file path separators in given path and returns the
 * position of the last one or std::string::npos if no file path separator was
 * found.
 */
static
std::string::size_type findLastPathSeparator(std::string const& path)
{
    return path.find_last_of("/\\");
}

/** Finds the position of last dot.
 * This could be used to extract file extension.
 */
static
std::string::size_type findLastDot(std::string const& path)
{
    return path.find_last_of('.');
}
} // end namespace

std::string dropFileExtension(std::string const& filename)
{
    // Look for dots in filename.
    auto const p = findLastDot(filename);
    if (p == std::string::npos)
        return filename;

    // Check position of the last path separator.
    auto const s = findLastPathSeparator(filename);
    if (s != std::string::npos && p < s)
        return filename;

    return filename.substr(0, p);
}

std::string extractBaseName(std::string const& pathname)
{
    auto const p = findLastPathSeparator(pathname);
    if (p == std::string::npos)
        return pathname;
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
    auto const p = findLastDot(str);
    if (p == std::string::npos)
        return std::string();
    return str.substr(p + 1);
}

bool hasFileExtension(std::string const& extension, std::string const& filename)
{
    return boost::iequals(extension, getFileExtension(filename));
}

static const char pathSeparator =
#ifdef _WIN32
                            '\\';
#else
                            '/';
#endif

std::string copyPathToFileName(const std::string &file_name,
                               const std::string &source)
{
    // check if file_name already contains a full path
    auto const pos = findLastPathSeparator(file_name);
    if (pos != std::string::npos)
        return file_name;

    if (source.empty())
        return file_name;
    if (source.back() != pathSeparator)
        return BaseLib::extractPath(source + pathSeparator).append(file_name);
    return BaseLib::extractPath(source).append(file_name);
}

std::string extractPath(std::string const& pathname)
{
    auto const pos = findLastPathSeparator(pathname);
    if (pos == std::string::npos)
        return "";
    return pathname.substr(0, pos + 1);
}

std::string appendPathSeparator(std::string const& path)
{
    if(findLastPathSeparator(path) == path.length() - 1)
        return path;
    return path + pathSeparator;
}

std::string joinPaths(std::string const& pathA, std::string const& pathB)
{
    if (pathA.empty())
        return pathB;

    if (pathB.empty())
        return pathA;

    if (pathB.front() == pathSeparator) {
        auto const tmpB = pathB.substr(1);
        return appendPathSeparator(pathA) + tmpB;
    } else {
        return appendPathSeparator(pathA) + pathB;
    }
}

} // end namespace BaseLib
