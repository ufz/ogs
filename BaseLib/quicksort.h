/**
 * \brief  Definition of the quicksort function.
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
#include <cassert>
#include <cstddef>
#include <iterator>
#include <vector>

namespace BaseLib
{

/// @pre {first1 <= last1 and the second iterator can be incremented
/// distance(first1, last1) times}
template <typename It1, typename It2>
void quicksort(It1 first1, It1 last1, It2 first2)
{
    using T1 = typename std::iterator_traits<It1>::value_type;
    using T2 = typename std::iterator_traits<It2>::value_type;

    std::vector<std::pair<T1, T2>> data;
    data.reserve(std::distance(first1, last1));
    std::transform(
        first1, last1, first2, std::back_inserter(data),
        [](T1 const& t1, T2 const& t2) { return std::make_pair(t1, t2); });

    // Sort data using first element of the pair.
    std::sort(begin(data), end(data),
              [](std::pair<T1, T2> const& a, std::pair<T1, T2> const& b) {
                  return (a.first < b.first);
              });

    // Unzip sorted data.
    for (auto const& pair : data)
    {
        *first1 = pair.first;
        *first2 = pair.second;
        ++first1;
        ++first2;
    }
}

/// @pre {end<=array.size() and perm.size()==array.size()}
template <typename T1, typename T2 = std::size_t>
void quicksort(T1* array, std::size_t beg, std::size_t end, T2* perm)
{
    assert(beg <= end);

    quicksort(array + beg, array + end, perm + beg);
}

template <typename T1, typename T2 = std::size_t>
void quicksort(std::vector<T1>& array, std::size_t beg, std::size_t end, std::vector<T2>& perm)
{
    assert (beg<=end);
    assert (end<=array.size());
    assert (perm.size()==array.size());

    quicksort(array.data(), beg, end, perm.data());
}

template <typename T1, typename T2 = std::size_t>
void quicksort(std::vector<T1*>& array, std::size_t beg, std::size_t end, std::vector<T2>& perm)
{
    assert (beg<=end);
    assert (end<=array.size());
    assert (perm.size()==array.size());

    // Zip input arrays.
    std::vector<std::pair<T1*, T2>> data;
    data.reserve(end-beg);
    std::transform(array.begin()+beg, array.begin()+end, perm.begin()+beg,
        std::back_inserter(data),
        [](T1* const& t1, T2 const& t2)
        {
            return std::make_pair(t1, t2);
        });

    // Sort data using first element of the pair.
    std::sort(data.begin(), data.end(),
        [](std::pair<T1*, T2> const& a, std::pair<T1*, T2> const& b)
        {
            return (*a.first < *b.first);
        });

    // Unzip sorted data.
    for (std::size_t i = 0; i < data.size(); i++)
    {
        array[beg+i] = data[i].first;
        perm[beg+i] = data[i].second;
    }
}

} // end namespace BaseLib
