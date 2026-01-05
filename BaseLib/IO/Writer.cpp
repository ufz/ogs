// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    out.precision(std::numeric_limits<double>::max_digits10);
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

int writeStringToFile(std::string_view content,
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
