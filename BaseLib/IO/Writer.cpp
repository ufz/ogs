/**
 * \file
 * \author Lars Bilke
 * \date   2012-02-13
 * \brief  Implementation of the Writer class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Writer.h"

#include <fstream>
#include <limits>

#include "BaseLib/Logging.h"

namespace BaseLib
{
namespace IO
{

Writer::Writer()
{
    out_.precision(std::numeric_limits<double>::digits10);
}

std::string Writer::writeToString()
{
    // Empty stream and clear error states.
    out_.str("");
    out_.clear();

    if (this->write())
    {
        return out_.str();
    }

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
            ERR("Could not open file '{:s}'!", filename);
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
    out_.precision(precision);
}
}  // namespace IO
}  // namespace BaseLib
