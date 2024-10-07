/**
 * \file
 * \author Lars Bilke
 * \date   Apr. 2010
 * \brief Filename manipulation routines implementation.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FileTools.h"

#include <spdlog/fmt/bundled/core.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/endian/conversion.hpp>
#include <filesystem>
#include <fstream>
#include <typeindex>
#include <unordered_map>

#include "BaseLib/Logging.h"
#include "Error.h"

namespace
{
/// The directory where the prj file resides.
std::string project_directory;

/// Whether the project directory has already been set.
bool project_directory_is_set = false;
}  // anonymous namespace

namespace BaseLib
{
bool isProjectDirectorySet()
{
    return project_directory_is_set;
}

/**
 * Returns true if given file exists.
 */
bool IsFileExisting(const std::string& strFilename)
{
    return std::filesystem::exists(std::filesystem::path(strFilename));
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
bool substituteKeyword(std::string& result,
                       std::string const& parenthesized_string,
                       std::string::size_type const begin,
                       std::string::size_type const end,
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
    // see https://fmt.dev/latest/syntax/#format-specification-mini-language
    if (b == 'e' || b == 'E' || b == 'f' || b == 'F' || b == 'g' || b == 'G')
    {
        type_specification[std::type_index(typeid(double))] = b;
        precision_specification.pop_back();
    }

    std::string const generated_fmt_string =
        "{" + precision_specification +
        type_specification[std::type_index(typeid(data))] + "}";
    result.replace(
        begin, end - begin + 1,
        fmt::vformat(generated_fmt_string, fmt::make_format_args(data)));

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
    auto const filename_path = std::filesystem::path(filename);
    return (filename_path.parent_path() / filename_path.stem()).string();
}

std::string extractBaseName(std::string const& pathname)
{
    return std::filesystem::path(pathname).filename().string();
}

std::string extractBaseNameWithoutExtension(std::string const& pathname)
{
    std::string basename = extractBaseName(pathname);
    return dropFileExtension(basename);
}

std::string getFileExtension(const std::string& path)
{
    return std::filesystem::path(path).extension().string();
}

bool hasFileExtension(std::string const& extension, std::string const& filename)
{
    return boost::iequals(extension, getFileExtension(filename));
}

std::string extractPath(std::string const& pathname)
{
    return std::filesystem::path(pathname).parent_path().string();
}

std::string joinPaths(std::string const& pathA, std::string const& pathB)
{
    return (std::filesystem::path(pathA) /= std::filesystem::path(pathB))
        .string();
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
    // TODO Remove these global vars. They are a possible source of errors when
    // invoking OGS from Python multiple times within a single session.
    project_directory = dir;
    project_directory_is_set = true;
}

void unsetProjectDirectory()
{
    project_directory.clear();
    project_directory_is_set = false;
}

void removeFile(std::string const& filename)
{
    bool const success =
        std::filesystem::remove(std::filesystem::path(filename));
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

bool createOutputDirectory(std::string const& dir)
{
    if (dir.empty())
    {
        return false;
    }

    std::error_code mkdir_err;
    if (std::filesystem::create_directories(dir, mkdir_err))
    {
        INFO("Output directory {:s} created.", dir);
    }
    else if (mkdir_err.value() != 0)
    {
        WARN("Could not create output directory {:s}. Error code {:d}, {:s}",
             dir, mkdir_err.value(), mkdir_err.message());
        return false;
    }
    return true;
}

template <typename T>
T readBinaryValue(std::istream& in)
{
    T v;
    in.read(reinterpret_cast<char*>(&v), sizeof(T));
    return v;
}

// explicit template instantiation
template float readBinaryValue<float>(std::istream&);
template double readBinaryValue<double>(std::istream&);

template <typename T>
std::vector<T> readBinaryArray(std::string const& filename, std::size_t const n)
{
    std::ifstream in(filename.c_str());
    if (!in)
    {
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

// explicit template instantiation
template std::vector<float> readBinaryArray<float>(std::string const&,
                                                   std::size_t const);

template <typename T>
std::vector<T> readBinaryVector(std::string const& filename)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file)
    {
        OGS_FATAL("File {:s} for curve definition not found", filename);
    }

    // Determine file size
    file.seekg(0, std::ios::end);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    size_t num_elements = size / sizeof(T);

    // Initialize vector with the right size
    std::vector<T> result(num_elements);

    // Read data directly into the vector
    if (!file.read(reinterpret_cast<char*>(result.data()), size))
    {
        OGS_FATAL("Could not read data from file {:s}.", filename);
    }

    if constexpr (std::endian::native != std::endian::little)
    {
        boost::endian::endian_reverse_inplace(result);
    }

    std::streamsize bytes_read = file.gcount();

    if (bytes_read != size)
    {
        OGS_FATAL(
            "Incomplete read: Expected {:d} bytes, but read {:d} bytes "
            "from "
            "file {:s}.",
            size, bytes_read, filename);
    }
    return result;
}

// explicit template instantiation
template std::vector<float> readBinaryVector<float>(std::string const&);
template std::vector<double> readBinaryVector<double>(std::string const&);

template <typename T>
void writeValueBinary(std::ostream& out, T const& val)
{
    out.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

// explicit template instantiation
template void writeValueBinary<std::size_t>(std::ostream&, std::size_t const&);

}  // end namespace BaseLib
