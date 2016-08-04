/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MAKEVECTORUNIQUE_H
#define MAKEVECTORUNIQUE_H

#include <algorithm>
#include <vector>

namespace BaseLib
{

/// Make the entries of the std::vector \c v unique. The remaining entries will
/// be sorted.
template <typename T>
void makeVectorUnique(std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    auto it = std::unique(v.begin(), v.end());
    v.erase(it, v.end());
}

} // end namespace BaseLib
#endif
