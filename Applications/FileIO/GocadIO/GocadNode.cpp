/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
