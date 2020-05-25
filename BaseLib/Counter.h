/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cstddef>

namespace BaseLib
{
template <typename X>
struct Counter
{
    Counter()
    {
        counter_value_++;
    }
    static std::size_t counter_value_;
};

template <typename X> std::size_t Counter<X>::counter_value_(0);

} // end namespace BaseLib
