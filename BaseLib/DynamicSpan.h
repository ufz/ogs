/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

namespace BaseLib
{
template <typename T>
struct DynamicSpan
{
    T* data;
    std::size_t size;
};
}  // namespace BaseLib
