// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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

Layer parseLayer(std::string const& line, std::vector<Region> const& regions);

}  // end namespace Gocad
}  // end namespace FileIO
