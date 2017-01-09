/**
 * \brief  Definition of the quicksort function.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef QUICKSORT_H_
#define QUICKSORT_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <vector>

namespace BaseLib
{

/// @pre {end<=array.size() and perm.size()==array.size()}
template <typename T1, typename T2 = std::size_t>
void quicksort(T1* array, std::size_t beg, std::size_t end, T2* perm)
{
    assert (beg <= end);

    // Zip input arrays.
    std::vector<std::pair<T1, T2>> data;
    data.reserve(end-beg);
    std::transform(array+beg, array+end, perm+beg,
        std::back_inserter(data),
        [](T1 const& t1, T2 const& t2)
        {
            return std::make_pair(t1, t2);
        });

    // Sort data using first element of the pair.
    std::sort(data.begin(), data.end(),
        [](std::pair<T1, T2> const& a, std::pair<T1, T2> const& b)
        {
            return (a.first < b.first);
        });

    // Unzip sorted data.
    for (std::size_t i = 0; i < data.size(); i++)
    {
        array[beg+i] = data[i].first;
        perm[beg+i] = data[i].second;
    }
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

#endif /* QUICKSORT_H_ */
