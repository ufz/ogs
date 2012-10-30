/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file FileTools.cpp
 *
 * Created on 2010-04-26 by Lars Bilke
 *
 */

#include "FileTools.h"

#include <sys/stat.h>

namespace BaseLib
{

/**
 * Returns true if given file exists. From http://www.techbytes.ca/techbyte103.html
 */
bool IsFileExisting(const std::string &strFilename)
{
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;

	// Attempt to get the file attributes
	intStat = stat(strFilename.c_str(),&stFileInfo);

	if(intStat == 0)
	{
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	}
	else
	{
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;
	}

	return(blnReturn);
}

/**
 * \brief truncate a file
 */
void truncateFile( std::string const& filename)
{
    std::ofstream ofs(filename.c_str(), std::ios_base::trunc);
    ofs.close();
}

size_t findLastPathSeparator(std::string const& path)
{
	return path.find_last_of("/\\");
}

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

std::string getSuffixFromPath(const std::string &str)
{
	std::string::size_type beg(str.find_last_of('.'));
	return str.substr(beg+1, str.length()-beg-1);
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
			std::string path = BaseLib::extractPath(source);
			return path.append(file_name);
		}
		else return std::string(file_name);
	}
	else return std::string(file_name);
}

std::string extractPath(std::string const& pathname)
{
	// extract path for reading external files
	size_t pos(pathname.rfind("/")); // linux, mac delimiter
	if (pos == std::string::npos) {
		pos = pathname.rfind("\\"); // windows delimiter
		if (pos == std::string::npos)
			pos = 0;
	}
	return pathname.substr(0, pos==0 ? pos : pos + 1);
}

} // end namespace BaseLib

