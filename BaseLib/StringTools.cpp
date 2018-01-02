/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-16
 * \brief  Implementation of string helper functions.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "StringTools.h"

#include <algorithm>
#include <cctype>
#include <cstdarg>
#include <cstdio>
#include <iomanip>

#include <logog/include/logog.hpp>
#include <boost/algorithm/string/replace.hpp>

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

std::list<std::string> splitString(const std::string &str, char delim)
{
    std::list<std::string> strList;
    std::stringstream ss(str);
    std::string item;
    while(getline(ss, item, delim))
        strList.push_back(item);
    return strList;
}

std::string replaceString(const std::string &searchString,
                          const std::string &replaceString,
                          std::string stringToReplace)
{
    boost::replace_all(stringToReplace, searchString, replaceString);
    return stringToReplace;
}

void trim(std::string &str, char ch)
{
    std::string::size_type pos = str.find_last_not_of(ch);
    if(pos != std::string::npos)
    {
        str.erase(pos + 1);
        pos = str.find_first_not_of(ch);
        if(pos != std::string::npos)
            str.erase(0, pos);
    }
    else
        str.erase(str.begin(), str.end());
}

void simplify(std::string &str)
{
    trim (str);
    str.erase(
        std::unique(str.begin(), str.end(), [](char a, char b) { return a == ' ' && b == ' '; }),
        str.end()
    );
}

std::string padLeft(std::string const& str, int maxlen, char ch)
{
    std::stringstream ss(str);
    ss << std::right << std::setw(maxlen) << std::setfill(ch) << str;
    return ss.str();
}

std::string const& tostring(std::string const& value)
{
    return value;
}

std::string format(const char* format_str, ... )
{
    va_list args;
    va_start(args, format_str);
    // get the number of chars to write
    va_list args_tmp;
    va_copy(args_tmp, args);
    int char_length = std::vsnprintf(nullptr, 0, format_str, args_tmp);
    va_end(args_tmp);
    // allocate buffer and store formatted output there
    std::vector<char> buffer(char_length + 1); // note +1 for null terminator
    vsnprintf(buffer.data(), buffer.size(), format_str, args);
    va_end(args);

    return std::string(buffer.data());
}

} // end namespace BaseLib
