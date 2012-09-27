/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file StringTools.h
 *
 * Created on 2010-06-16 by Thomas Fischer
 */

#ifndef STRINGTOOLS_H
#define STRINGTOOLS_H

#include <string>
#include <list>
#include <sstream>
#include <fstream>
#include <iostream>
#include <ctype.h>


/**
 *   Splits a string into a list of strings.
 *  \param str String to be splitted
 *  \param delim Character indicating that the string should be splitted
 *  \return
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

namespace BaseLib {

/**
 * Extract the filename from a path
 */
std::string getFileNameFromPath(const std::string &str, bool with_extension = false);

/**
 * Extract the file type / suffix from a path
 */
std::string getSuffixFromPath(const std::string &str);


/**
 * Checks if file_name already contains a qualified path and if not copies the path from source.
 */
std::string copyPathToFileName(const std::string &file_name, const std::string &source);

/**
 * extracts the path of a fully qualified path name of the file
 * @param fname [input] the fully qualified path name of the file
 * @param path [output] the path of the fully qualified path name of the file
 */
void extractPath (std::string const& fname, std::string& path);


} // end namespace BaseLib

#ifdef MSVC
void correctScientificNotation(std::string filename, std::size_t precision = 0);
#endif

#endif //STRINGTOOLS_H
