// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Region.h"

#include <iterator>
#include <sstream>

#include "BaseLib/Logging.h"

namespace FileIO
{
namespace Gocad
{
std::ostream& operator<<(std::ostream& os, Region const& r)
{
    return os << "(" << r.name << "|" << r.bit << ")";
}

Region parseRegion(std::string const& line)
{
    std::istringstream iss(line);
    std::istream_iterator<std::string> it(iss);
    // Check first word is REGION or MODEL_REGION.
    if (*it != std::string("REGION") && *it != std::string("MODEL_REGION"))
    {
        ERR("Expected REGION or MODEL_REGION keyword but '{:s}' found.\n",
            it->c_str());
        throw std::runtime_error(
            "In parseRegion() expected REGION or MODEL_REGION keyword not "
            "found.\n");
    }
    ++it;

    Region r;
    r.name = *it;
    ++it;
    r.bit = atoi(it->c_str());

    return r;
}

}  // end namespace Gocad
}  // end namespace FileIO
