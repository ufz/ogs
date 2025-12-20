// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "readStringListFromFile.h"

#include <fstream>

#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"

namespace BaseLib
{
namespace IO
{
std::vector<std::string> readStringListFromFile(std::string const& filename)
{
    std::vector<std::string> string_list;
    std::ifstream in(filename);
    if (!in)
    {
        ERR("Could not open file {:s}.", filename);
        return string_list;
    }
    std::string line;
    while (std::getline(in, line))
    {
        trim(line);
        if (line.empty())
        {
            continue;
        }
        string_list.push_back(line);
    }
    return string_list;
}
}  // namespace IO
}  // namespace BaseLib
