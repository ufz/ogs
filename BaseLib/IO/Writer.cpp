/**
 * \file
 * \author Lars Bilke
 * \date   2012-02-13
 * \brief  Implementation of the Writer class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
    out.precision(std::numeric_limits<double>::digits10);
}

std::string Writer::writeToString()
{
    // Empty stream and clear error states.
    out.str("");
    out.clear();

    if (this->write())
    {
        return out.str();
    }

    return std::string("");
}

int writeStringToFile(std::string content,
                      std::filesystem::path const& file_path)
{
    if (content.empty())
    {
        return 0;
    }
    std::ofstream fileStream;
    fileStream.open(file_path.c_str());

    // check file stream
    if (!fileStream)
    {
        ERR("Could not open file '{:s}'!", file_path.string());
        return 0;
    }

    fileStream << content;

    fileStream.close();
    return 1;
}

}  // namespace IO
}  // namespace BaseLib
