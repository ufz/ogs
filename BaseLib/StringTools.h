/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-16
 * \brief  Definition of string helper functions.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef STRINGTOOLS_H
#define STRINGTOOLS_H


#include <string>
#include <list>
#include <sstream>
#include <vector>

namespace BaseLib {

/**
 *  Splits a string into a vector of strings. This method only works for string seperation 
 *  recognised by the std::stringstream iterator such as ' ' or '\t'.
 *  \param str String to be splitted
 *  \return Vector of strings
 */
std::vector<std::string> splitString(std::string const& str);

/**
 *   Splits a string into a list of strings.
 *  \param str String to be splitted
 *  \param delim Character indicating that the string should be splitted
 *  \return List of strings
 */
std::list<std::string> splitString(const std::string &str, char delim);

/**
 *   Replaces a substring with another in a string
 *  \param searchString Search for this string
 *  \param replaceString Replace with this string
 *  \param stringToReplace Search and replace in this string
 *  \return The modified string
 */
std::string replaceString(const std::string &searchString, const std::string &replaceString, std::string stringToReplace);

/**
 *   Converts a string into a number (double, float, int, ...)
 *  Example: std::size_t number (str2number<std::size_t> (str));
 *  \param str string to be converted
 *  \return the number
 */
template<typename T> T str2number (const std::string &str)
{
    std::stringstream strs (str, std::stringstream::in | std::stringstream::out);
    T v;
    strs >> v;
    return v;
}

/**
 * Strip whitespace (or other characters) from the beginning and end of a string.
 * Equivalent functionality to Qt::QString::trim().
 */
void trim(std::string &str, char ch=' ');

/**
 * Removes multiple whitespaces (or other characters) from within a string.
 * Equivalent functionality to Qt::QString::simplify().
 */
void simplify(std::string &str);

/**
 * Returns the string which is right aligned with padding on the left.
 */
std::string padLeft(std::string const& str, int maxlen, char ch=' ');


//! Method for handling conversion to string uniformly across all types and std::string; see std::string overload below.
template<typename T> std::string tostring(T const& value)
{
    return std::to_string(value);
}
//! \overload
std::string const& tostring(std::string const& value);

//! returns printf-like formatted string
std::string format(const char* format_string, ... );

} // end namespace BaseLib

#ifdef MSVC
void correctScientificNotation(std::string filename, std::size_t precision = 0);
#endif

#endif //STRINGTOOLS_H
