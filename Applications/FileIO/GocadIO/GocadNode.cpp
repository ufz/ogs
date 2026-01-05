// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "GocadNode.h"

namespace FileIO
{
namespace Gocad
{
bool operator<=(GocadNode const& n0, GocadNode const& n1)
{
    for (std::size_t k(0); k < 3; k++)
    {
        if (n0[0] > n1[0])
        {
            return false;
        }
        if (n0[0] < n1[0])
        {
            return true;
        }
        // => n0[k] == n1[k]
    }

    return n0.getLayerTransitionIndex() <= n1.getLayerTransitionIndex();
}

}  // end namespace Gocad
}  // end namespace FileIO
