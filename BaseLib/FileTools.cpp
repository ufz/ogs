/**
 * \file
 * \author Lars Bilke
 * \date   Apr. 2010
 * \brief Filename manipulation routines implementation.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FileTools.h"

#include <boost/algorithm/string.hpp>
#include <typeindex>
#include <unordered_map>

#include "Error.h"
#include "filesystem.h"

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
 * Returns true if given file exists.
 */
bool IsFileExisting(const std::string& strFilename)
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

std::string containsKeyword(std::string const& str, std::string const& keyword)
{
    auto const position = str.find(keyword);
    if (position != std::string::npos)
    {
        return str.substr(0, position);
    }
    return "";
}

template <typename T>
bool substituteKeyword(std::string& result, std::string& parenthesized_string,
                       std::string::size_type begin, std::string::size_type end,
                       std::string const& keyword, T& data)
{
    std::string precision_specification =
        containsKeyword(parenthesized_string, keyword);

    if (precision_specification.empty())
    {
        return false;
    }

    std::unordered_map<std::type_index, char> type_specification;
    type_specification[std::type_index(typeid(int))] = 'd';
    type_specification[std::type_index(typeid(double))] = 'f';  // default
    type_specification[std::type_index(typeid(std::string))] = 's';

    auto const& b = precision_specification.back();
    // see https://fmt.dev/latest/syntax.html#formatspec
    if (b == 'e' || b == 'E' || b == 'f' || b == 'F' || b == 'g' || b == 'G')
    {
        type_specification[std::type_index(typeid(double))] = b;
        precision_specification = precision_specification.substr(
            0, precision_specification.length() - 1);
    }

    std::string const generated_fmt_string =
        "{" + precision_specification +
        type_specification[std::type_index(typeid(data))] + "}";
    result = result.substr(0, begin) + fmt::format(generated_fmt_string, data) +
             result.substr(end + 1, result.length() - (end + 1));
    return true;
}

std::string constructFormattedFileName(std::string const& format_specification,
                                       std::string const& mesh_name,
                                       int const timestep,
                                       double const t,
                                       int const iteration)
{
    char const open_char = '{';
    char const close_char = '}';
    std::string::size_type begin = 0;
    std::string::size_type end = std::string::npos;
    std::string result = format_specification;

    while (begin != std::string::npos)
    {
        auto length_before_substitution = result.length();
        // find next parenthesized string
        std::string str = "";
        std::tie(str, begin, end) =
            getParenthesizedString(result, open_char, close_char, begin);
        if (!substituteKeyword(result, str, begin, end, "timestep", timestep) &&
            !substituteKeyword(result, str, begin, end, "time", t) &&
            !substituteKeyword(result, str, begin, end, "iteration", iteration))
        {
            substituteKeyword(result, str, begin, end, "meshname", mesh_name);
        }
        begin = end - (length_before_substitution - result.length());
    }

    return result;
}

double swapEndianness(double const& v)
{
    union
    {
        double v;
        char c[sizeof(double)];
    } a{}, b{};

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

std::string getFileExtension(const std::string& path)
{
    return fs::path(path).extension().string();
}

bool hasFileExtension(std::string const& extension, std::string const& filename)
{
    return boost::iequals(extension, getFileExtension(filename));
}

std::string copyPathToFileName(const std::string& file_name,
                               const std::string& source)
{
    auto filePath = fs::path(file_name);
    if (filePath.has_parent_path())
    {
        return filePath.string();
    }
    return (fs::path(source) /= filePath).string();
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

void removeFile(std::string const& filename)
{
    bool const success = fs::remove(fs::path(filename));
    if (success)
    {
        DBUG("Removed '{:s}'", filename);
    }
}

void removeFiles(std::vector<std::string> const& files)
{
    for (auto const& file : files)
    {
        removeFile(file);
    }
}
}  // end namespace BaseLib
