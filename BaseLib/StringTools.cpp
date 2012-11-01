/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file StringTools.cpp
 *
 * Created on 2010-06-16 by Thomas Fischer
 */

#include "StringTools.h"

#include <algorithm>
#include <cctype>

namespace BaseLib {

std::list<std::string> splitString(const std::string &str, char delim)
{
	std::list<std::string> strList;
	std::stringstream ss(str);
    std::string item;
    while(getline(ss, item, delim)) {
        strList.push_back(item);
    }
    return strList;
}

std::string replaceString(const std::string &searchString, const std::string &replaceString, std::string stringToReplace)
 {
	std::string::size_type pos = stringToReplace.find(searchString, 0);
	int intLengthSearch = searchString.length();

	while (std::string::npos != pos) {
		stringToReplace.replace(pos, intLengthSearch, replaceString);
		pos = stringToReplace.find(searchString, 0);
	}
	return stringToReplace;
}

void trim(std::string &str, char ch)
{
  std::string::size_type pos = str.find_last_not_of(ch);
  if(pos != std::string::npos)
  {
    str.erase(pos + 1);
    pos = str.find_first_not_of(ch);
    if(pos != std::string::npos) str.erase(0, pos);
  }
  else str.erase(str.begin(), str.end());
}

std::string stringToUpper(std::string const& str)
{
    std::string s = str;
	std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) std::toupper);
    return s;
}

} // end namespace BaseLib


#ifdef MSVC
void correctScientificNotation(std::string filename, size_t precision)
{
	std::ifstream stream;
	std::ofstream outputStream;

	stream.open(filename.c_str());
	std::string tmpFilename = filename + ".tmp";
	outputStream.open(tmpFilename.c_str());

	if (!stream)
	{
		std::cout << "correctScientificNotation: fstream is not open" << std::endl;
		return;
	}

	std::string line;

	// Iterate over lines in stream
	while (getline(stream, line))
	{
		std::string word;
		std::istringstream iss(line);
		// Iterate over all words in line
		while (iss >> word)
		{
			// Search for e+0
			std::size_t exponentPosition = word.find("e+0", precision);
			if (exponentPosition == std::string::npos)
				// If not found search for e-0
				exponentPosition = word.find("e-0", precision);
			if (exponentPosition != std::string::npos)
			{
				std::size_t wordSize = word.size();
				std::size_t exponentSize = wordSize - exponentPosition;

				if(exponentSize > 4)
				{
					// Erase the leading zero considering trailing characters
					int i = wordSize - 1;
					while (!isdigit(word[i]))
						--i;

					size_t erasePos = wordSize - 3 - (wordSize - 1 - i);
					std::string eraseString = word.substr(erasePos, 1);
					if (eraseString.find("0") != std::string::npos)
						word.erase(erasePos, 1);
				}
			}

			outputStream << word << " ";
		}
		outputStream << std::endl;
	}

	stream.close();
	outputStream.close();

	remove(filename.c_str());
	rename(tmpFilename.c_str(), filename.c_str());
}
#endif


