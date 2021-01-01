/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <algorithm>
#include <string>
#include <vector>

#include "Region.h"

namespace FileIO
{
namespace Gocad
{
// Each model layer own several regions.
struct Layer final
{
    std::vector<Region> regions;

    bool hasRegion(Region const& r) const
    {
        return (std::find(regions.begin(), regions.end(), r) != regions.end());
    }
};

std::ostream& operator<<(std::ostream& os, Layer const& l);

Layer parseLayer(std::string const& line,
                 std::vector<Gocad::Region> const& regions);

}  // end namespace Gocad
}  // end namespace FileIO
