/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

namespace NumLib
{

template <typename IndexType>
struct IndexValueVector final
{
    std::vector<IndexType> ids;
    std::vector<double> values;
};

}
