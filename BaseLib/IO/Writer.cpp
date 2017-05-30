/**
 * \file
 * \author Lars Bilke
 * \date   2012-02-13
 * \brief  Implementation of the Writer class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Writer.h"

#include <fstream>
#include <limits>

#include <logog/include/logog.hpp>

namespace BaseLib
{
namespace IO
{

Writer::Writer()
{
    _out.precision(std::numeric_limits<double>::digits10);
}

std::string Writer::writeToString()
{
    // Empty stream and clear error states.
    _out.str("");
    _out.clear();

    if (this->write())
        return _out.str();

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
            ERR("Could not open file \"%s\"!", filename.c_str());
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

} // namespace IO
} // namespace BaseLib
