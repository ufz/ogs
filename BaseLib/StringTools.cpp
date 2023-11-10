/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-16
 * \brief  Implementation of string helper functions.
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "StringTools.h"

#include <algorithm>
#include <boost/algorithm/string/replace.hpp>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <iterator>
#include <random>

#include "Error.h"

namespace BaseLib
{
std::vector<std::string> splitString(std::string const& str)
{
    std::istringstream str_stream(str);
    std::vector<std::string> items;
    std::copy(std::istream_iterator<std::string>(str_stream),
              std::istream_iterator<std::string>(),
              std::back_inserter(items));
    return items;
}

std::list<std::string> splitString(const std::string& str, char delim)
{
    std::list<std::string> strList;
    std::stringstream ss(str);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        strList.push_back(item);
    }
    return strList;
}

std::string replaceString(const std::string& searchString,
                          const std::string& replaceString,
                          std::string stringToReplace)
{
    boost::replace_all(stringToReplace, searchString, replaceString);
    return stringToReplace;
}

void trim(std::string& str, char ch)
{
    std::string::size_type pos = str.find_last_not_of(ch);
    if (pos != std::string::npos)
    {
        str.erase(pos + 1);
        pos = str.find_first_not_of(ch);
        if (pos != std::string::npos)
        {
            str.erase(0, pos);
        }
    }
    else
    {
        str.erase(str.begin(), str.end());
    }
}

void simplify(std::string& str)
{
    trim(str);
    str.erase(std::unique(str.begin(),
                          str.end(),
                          [](char a, char b) { return a == ' ' && b == ' '; }),
              str.end());
}

std::string format(const char* format_str, ...)
{
    va_list args;
    va_start(args, format_str);
    // get the number of chars to write
    va_list args_tmp;
    va_copy(args_tmp, args);
    int char_length = std::vsnprintf(nullptr, 0, format_str, args_tmp);
    va_end(args_tmp);
    // allocate buffer and store formatted output there
    std::vector<char> buffer(char_length + 1);  // note +1 for null terminator
    vsnprintf(buffer.data(), buffer.size(), format_str, args);
    va_end(args);

    return std::string(buffer.data());
}

std::string randomString(std::size_t const length)
{
    static constexpr char charset[] =
        "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

    static const auto seed = static_cast<std::mt19937::result_type>(
        std::chrono::system_clock::now().time_since_epoch().count());
    static std::mt19937 generator{seed};
    static std::uniform_int_distribution<unsigned short> distribution(
        0, sizeof(charset) - 2);

    std::string s(length, 0);
    std::generate_n(begin(s), length,
                    [&]() { return charset[distribution(generator)]; });

    return s;
}

std::string getUniqueName(std::vector<std::string> const& existing_names,
                          std::string const& input_name)
{
    std::string result_name = input_name;
    std::size_t count = 1;
    while (std::find(existing_names.cbegin(),
                     existing_names.cend(),
                     result_name) != existing_names.end())
    {
        count++;
        result_name = input_name + "-" + std::to_string(count);
    }
    return result_name;
}

std::vector<int> splitMaterialIdString(std::string const& material_id_string)
{
    auto const material_ids_strings =
        BaseLib::splitString(material_id_string, ',');

    std::vector<int> material_ids;
    for (auto& mid_str : material_ids_strings)
    {
        std::size_t num_chars_processed = 0;
        int material_id;
        try
        {
            material_id = std::stoi(mid_str, &num_chars_processed);
        }
        catch (std::invalid_argument&)
        {
            OGS_FATAL(
                "Could not parse material ID from '{}' to a valid "
                "integer.",
                mid_str);
        }
        catch (std::out_of_range&)
        {
            OGS_FATAL(
                "Could not parse material ID from '{}'. The integer value "
                "of the given string exceeds the permitted range.",
                mid_str);
        }

        if (num_chars_processed != mid_str.size())
        {
            // Not the whole string has been parsed. Check the rest.
            if (auto const it = std::find_if_not(
                    begin(mid_str) + num_chars_processed, end(mid_str),
                    [](unsigned char const c) { return std::isspace(c); });
                it != end(mid_str))
            {
                OGS_FATAL(
                    "Could not parse material ID from '{}'. Please "
                    "separate multiple material IDs by comma only. "
                    "Invalid character: '{}' at position {}.",
                    mid_str, *it, distance(begin(mid_str), it));
            }
        }

        material_ids.push_back(material_id);
    };

    return material_ids;
}

}  // end namespace BaseLib
