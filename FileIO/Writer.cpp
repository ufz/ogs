/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file Writer.cpp
 *
 * Created on 2012-02-13 by Lars Bilke
 */

// ** INCLUDES **
#include "Writer.h"

#include <fstream>

namespace FileIO
{

Writer::Writer()
{
}

std::string Writer::writeToString()
{
	// Empty stream and clear error states.
	_out.str("");
	_out.clear();

	if (this->write(_out))
		return _out.str();
	else
		return std::string("");
}

int Writer::writeToFile(std::string const& filename)
{
	std::string file_content = this->writeToString();
	if (!file_content.empty())
	{
		std::ofstream fileStream;
		fileStream.open (filename.c_str());

		// check file stream
		if (!fileStream)
		{
			std::cerr << "Could not open file " << filename << " !" << std::endl;
			return 0;
		}

		fileStream << file_content;

		fileStream.close();
		return 1;
	}
	return 0;
}

void Writer::setPrecision(unsigned int precision)
{
	_out.precision(precision);
}

void Writer::setFormat(std::ios_base::fmtflags flags)
{
	_out.setf(flags);
}

} // namespace FileIO
