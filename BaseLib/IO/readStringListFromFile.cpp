/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/IO/readStringListFromFile.h"

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
