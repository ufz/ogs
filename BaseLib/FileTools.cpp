/**
 * \file
 * \author Lars Bilke
 * \date   Apr. 2010
 * \brief Filename manipulation routines implementation.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FileTools.h"
#include "Error.h"
#include "filesystem.h"

#include <boost/algorithm/string.hpp>

namespace
{
/// The directory where the prj file resides.
std::string project_directory;

/// Whether the project directory has already been set.
bool project_directory_is_set = false;
}  // anonymous namespace

namespace BaseLib
{
/**
 * Returns true if given file exists. From http://www.techbytes.ca/techbyte103.html
 */
bool IsFileExisting(const std::string &strFilename)
{
    return fs::exists(fs::path(strFilename));
}

std::tuple<std::string, std::string::size_type, std::string::size_type>
getParenthesizedString(std::string const& in,
                       char const open_char,
                       char const close_char,
                       std::string::size_type pos)
{
    auto const pos_curly_brace_open = in.find_first_of(open_char, pos);
    if (pos_curly_brace_open == std::string::npos)
    {
        return std::make_tuple("", std::string::npos, std::string::npos);
    }
    auto const pos_curly_brace_close =
        in.find_first_of(close_char, pos_curly_brace_open);
    if (pos_curly_brace_close == std::string::npos)
    {
        return std::make_tuple("", std::string::npos, std::string::npos);
    }
    return std::make_tuple(
        in.substr(pos_curly_brace_open + 1,
                  pos_curly_brace_close - (pos_curly_brace_open + 1)),
        pos_curly_brace_open, pos_curly_brace_close);
}
std::string constructFileName(std::string const& prefix,
                              int const process_id,
                              int const timestep,
                              double const t)
{
    return prefix + "_pcs_" + std::to_string(process_id) + "_ts_" +
           std::to_string(timestep) + "_t_" + std::to_string(t);
}

double swapEndianness(double const& v)
{
    union
    {
        double v;
        char c[sizeof(double)];
    } a {}, b {};

    a.v = v;
    for (unsigned short i = 0; i < sizeof(double) / 2; i++)
    {
        b.c[i] = a.c[sizeof(double) / 2 - i - 1];
    }

    for (unsigned short i = sizeof(double) / 2; i < sizeof(double); i++)
    {
        b.c[i] = a.c[sizeof(double) + sizeof(double) / 2 - i - 1];
    }

    return b.v;
}

std::string dropFileExtension(std::string const& filename)
{
    auto const filename_path = fs::path(filename);
    return (filename_path.parent_path() / filename_path.stem()).string();
}

std::string extractBaseName(std::string const& pathname)
{
    return fs::path(pathname).filename().string();
}

std::string extractBaseNameWithoutExtension(std::string const& pathname)
{
    std::string basename = extractBaseName(pathname);
    return dropFileExtension(basename);
}

std::string getFileExtension(const std::string &path)
{
    return fs::path(path).extension().string();
}

bool hasFileExtension(std::string const& extension, std::string const& filename)
{
    return boost::iequals(extension, getFileExtension(filename));
}

std::string copyPathToFileName(const std::string &file_name,
                               const std::string &source)
{
    auto filePath = fs::path(file_name);
    if(filePath.has_parent_path())
    {
        return filePath.string();
    }
    else
    {
        return (fs::path(source) /= filePath).string();
    }

}

std::string extractPath(std::string const& pathname)
{
    return fs::path(pathname).parent_path().string();
}

std::string joinPaths(std::string const& pathA, std::string const& pathB)
{
    return (fs::path(pathA) /= fs::path(pathB)).string();
}

std::string const& getProjectDirectory()
{
    if (!project_directory_is_set)
    {
        OGS_FATAL("The project directory has not yet been set.");
    }
    return project_directory;
}

void setProjectDirectory(std::string const& dir)
{
    if (project_directory_is_set)
    {
        OGS_FATAL("The project directory has already been set.");
    }
    project_directory = dir;
    project_directory_is_set = true;
}

void removeFiles(std::vector<std::string> const& files)
{
    for (auto const& file : files)
    {
        bool const success = fs::remove(fs::path(file));
        if (success)
        {
            DBUG("Removed '{:s}'", file);
        }
    }
}
} // end namespace BaseLib
