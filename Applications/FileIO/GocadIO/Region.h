/**
 *
 * @copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#pragma once

#include <string>
#include <vector>

namespace FileIO
{
namespace Gocad
{
struct Region final
{
    std::string name;
    unsigned bit{};

    bool operator==(Region const& r) const { return bit == r.bit; }
};

std::ostream& operator<<(std::ostream& os, Region const& r);

Region parseRegion(std::string const& line);

}  // end namespace Gocad
}  // end namespace FileIO
