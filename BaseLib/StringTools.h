/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-16
 * \brief  Definition of string helper functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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
#include <fstream>
#include <iostream>
#include <ctype.h>

namespace BaseLib {

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
 *   Converts a number (double, float, int, ...) into a string
 *  \param d The number to be converted
 *  \return The number as string
 */
template<typename T> std::string number2str(T d)
{
	std::stringstream out;
	out << d;
	return out.str();
}

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
 */
void trim(std::string &str, char ch=' ');

/**
 * Returns same string with all characters in upper case.
 *
 * This uses std::toupper() function, and does not care about unicode.
 */
std::string stringToUpper(std::string const& str);

} // end namespace BaseLib

#ifdef MSVC
void correctScientificNotation(std::string filename, std::size_t precision = 0);
#endif

#endif //STRINGTOOLS_H
