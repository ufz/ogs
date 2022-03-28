/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <cstddef>

namespace BaseLib
{
template <typename T>
struct DynamicSpan
{
    template <typename Iterator>
    DynamicSpan(Iterator data_, std::size_t size)
        : data{&*data_ /* manually convert iterator to pointer */}, size_{size}
    {
    }

    T* begin() const { return data; }
    T* end() const { return data + size_; }
    std::size_t size() const { return size_; }
    T& operator[](std::size_t i) const { return data[i]; }

    T* data;

private:
    std::size_t size_;
};

template <typename T>
DynamicSpan(T*, std::size_t) -> DynamicSpan<T>;
}  // namespace BaseLib
