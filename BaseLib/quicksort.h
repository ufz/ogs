/**
 * \brief  Definition of the quicksort function.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef QUICKSORT_H_
#define QUICKSORT_H_

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <vector>

namespace BaseLib
{

template <typename T1, typename T2 = std::size_t>
void quicksort(std::vector<T1>& array, std::size_t beg, std::size_t end, std::vector<T2>& perm)
{
	// Zip input arrays.
	std::vector<std::pair<T1, T2>> data;
	data.reserve(end-beg);
	std::transform(array.begin()+beg, array.begin()+end, perm.begin()+beg,
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
	for (std::size_t i = beg; i < end; i++)
	{
		array[i] = data[i-beg].first;
		perm[i] = data[i-beg].second;
	}
}

template <typename T1, typename T2 = std::size_t>
void quicksort(std::vector<T1*>& array, std::size_t beg, std::size_t end, std::vector<T2>& perm)
{
	// Zip input arrays.
	std::vector<std::pair<T1*, T2>> data;
	data.reserve(end-beg);
	std::transform(array.begin()+beg, array.begin()+(end-beg), perm.begin()+beg,
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

template <typename T1, typename T2 = std::size_t>
void quicksort(T1* array, std::size_t beg, std::size_t end, T2* perm)
{
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

} // end namespace BaseLib

#endif /* QUICKSORT_H_ */
