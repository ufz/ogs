/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-16
 * \brief  Definition of string helper functions.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <list>
#include <sstream>
#include <string>
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
}  // end namespace BaseLib
