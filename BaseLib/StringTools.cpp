/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file StringTools.cpp
 *
 * Created on 2010-06-16 by Thomas Fischer
 */

#include "StringTools.h"

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

namespace BaseLib {

std::string getFileNameFromPath(const std::string &str, bool with_extension)
{
	std::string::size_type beg1 = str.find_last_of('/');
	std::string::size_type beg2 = str.find_last_of('\\');
	std::string::size_type beg;
	if (beg1 == std::string::npos && beg2 == std::string::npos) beg = -1;
	else if (beg1 == std::string::npos) beg = beg2;
	else if (beg2 == std::string::npos) beg = beg1;
	else beg = (beg1<beg2) ? beg2 : beg1;
	std::string file ( str.substr(beg+1) );
	if (with_extension) return file;
	// cut extension
	std::string::size_type end  = file.find_last_of('.');
	return file.substr(0,end);
}

std::string copyPathToFileName(const std::string &file_name, const std::string &source)
{
	// check if file_name already contains a full path
	size_t pos(file_name.rfind("/")); // linux, mac delimiter
	if (pos == std::string::npos)
	{
		pos = file_name.rfind("\\"); // windows delimiter
		if (pos == std::string::npos)
		{
			std::string path;
			BaseLib::extractPath(source, path);
			return path.append(file_name);
		}
		else return std::string(file_name);
	}
	else return std::string(file_name);
}


void extractPath (std::string const& fname, std::string& path)
{
	// extract path for reading external files
	size_t pos(fname.rfind("/")); // linux, mac delimiter
	if (pos == std::string::npos) {
		pos = fname.rfind("\\"); // windows delimiter
		if (pos == std::string::npos)
			pos = 0;
	}
	path = fname.substr(0, pos==0 ? pos : pos + 1);
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


