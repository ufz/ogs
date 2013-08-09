/**
 * \file
 * \author Lars Bilke
 * \date   2012-02-13
 * \brief  Implementation of the Writer class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "Writer.h"

#include <fstream>
#include <limits>

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
	_out.precision(std::numeric_limits<double>::digits10);

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
			std::cerr << "Could not open file " << filename << " !\n";
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
