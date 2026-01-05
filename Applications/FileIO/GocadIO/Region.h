// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
