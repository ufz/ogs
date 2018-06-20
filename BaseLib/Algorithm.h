/**
 * \file
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
#include "Error.h"

namespace BaseLib
{
template <typename InputIt, typename Predicate>
typename std::iterator_traits<InputIt>::reference findElement(
    InputIt begin, InputIt end, Predicate predicate,
    std::string const& error = "")
{
    auto it = std::find_if(begin, end, predicate);
    if (it == end)
    {
        OGS_FATAL("Element not found in the input range; %s", error.c_str());
    }
    return *it;
}

}  // end namespace BaseLib
