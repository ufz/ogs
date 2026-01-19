// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <list>
#include <optional>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/istream.hpp>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace BaseLib
{
/**
 *  Splits a string into a vector of strings. This method only works for string
 *  separation recognised by the std::stringstream iterator such as ' ' or
 * '\\t'.
 *  \param str String to be split
 *  \return Vector of strings
 */
std::vector<std::string> splitString(std::string const& str);

/**
 *   Splits a string into a list of strings.
 *  \param str String to be split
 *  \param delim Character indicating that the string should be split
 *  \return List of strings
 */
std::list<std::string> splitString(const std::string& str, char delim);

/**
 *   Replaces a substring with another in a string
 *  \param searchString Search for this string
 *  \param replaceString Replace with this string
 *  \param stringToReplace Search and replace in this string
 *  \return The modified string
 */
std::string replaceString(const std::string& searchString,
                          const std::string& replaceString,
                          std::string stringToReplace);

/**
 *   Converts a string into a number (double, float, int, ...)
 *  Example: std::size_t number (str2number<std::size_t> (str));
 *  \param str string to be converted
 *  \return the number
 */
template <typename T>
T str2number(const std::string& str)
{
    std::stringstream strs(str, std::stringstream::in | std::stringstream::out);
    T v;
    strs >> v;
    return v;
}

/**
 * Strip whitespace (or other characters) from the beginning and end of a
 * string. Equivalent functionality to Qt::QString::trim().
 */
void trim(std::string& str, char ch = ' ');

/**
 * Removes multiple whitespaces (or other characters) from within a string.
 * Equivalent functionality to Qt::QString::simplify().
 */
void simplify(std::string& str);

//! Returns a random string of the given length containing just a-z,A-Z,0-9
std::string randomString(std::size_t length);

//! Append '-' and a number such that the name is unique.
std::string getUniqueName(std::vector<std::string> const& existing_names,
                          std::string const& input_name);

/**
 *   Tries to parse a whitespace-separated list of values from a string.
 *   Returns std::nullopt if parsing fails and sets the bad token index.
 */
template <typename T>
std::optional<std::vector<T>> tryParseVector(std::string const& raw,
                                             std::size_t* bad_token_idx)
{
    std::istringstream iss{raw};

    // Create a range that reads T values from the stream
    auto values = ranges::istream_view<T>(iss);
    std::vector<T> out = ranges::to<std::vector>(values);

    // Check if we consumed the entire input
    if (!iss.eof())
    {
        if (bad_token_idx)
        {
            *bad_token_idx = out.size() + 1;
        }
        return std::nullopt;
    }
    return out;
}

}  // end namespace BaseLib
