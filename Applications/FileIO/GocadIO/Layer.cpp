/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Layer.h"

#include <algorithm>
#include <iterator>
#include <sstream>

#include "BaseLib/Logging.h"

namespace FileIO
{
namespace Gocad
{

std::ostream& operator<<(std::ostream& os, Layer const& l)
{
    std::copy(l.regions.begin(), l.regions.end(),
              std::ostream_iterator<Gocad::Region>(os, " "));
    return os;
}


Layer parseLayer(std::string const& line, std::vector<Gocad::Region> const& regions)
{
    std::istringstream iss(line);
    std::istream_iterator<std::string> it(iss);
    // Check first word is MODEL_LAYER.
    if (*it != std::string("MODEL_LAYER"))
    {
        ERR("Expected MODEL_LAYER keyword but '{:s}' found.\n", it->c_str());
        throw std::runtime_error(
            "In parseRegion() expected MODEL_LAYER keyword not found.\n");
    }
    ++it;

    Layer l;
    while (it != std::istream_iterator<std::string>() && *it != "END")
    {
        auto const& region_it =
            std::find_if(regions.begin(), regions.end(),
                         [&](Gocad::Region const& r) { return r.name == *it; });
        if (region_it != regions.end())
        {
            l.regions.push_back(*region_it);
        }
        ++it;
    }

    return l;
}

}  // end namespace Gocad
}  // end namespace FileIO
