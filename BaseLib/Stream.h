/**
 * \author Norihiro Watanabe
 * \date   2014-03-07
 * \brief  Helper tools for debugging
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <iterator>
#include <ostream>
#include <vector>

template<typename T>
std::ostream &operator <<(std::ostream &os, const std::vector<T> &v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    os << "\n";
    return os;
}
