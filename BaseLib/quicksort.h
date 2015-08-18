/**
 * \file
 * \author Thomas Fischer
 * \date   2010-05-26
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

// STL
#include <algorithm>
#include <cstddef>

namespace BaseLib
{

/**
 * Permutes the entries of a part of an array such that all entries that are smaller
 * than a certain value are at the beginning of the array and all entries that are
 * bigger are at the end of the array. This version of partition_ permutes a second
 * array second_array according to the sorting.
 * @param array array to sort
 * @param beg beginning index in array for sorting
 * @param end end-1 is the last index in array for sorting
 * @param second_array the second array is permuted according to the sort process of array
 */
template <typename T1, typename T2>
std::size_t partition_(T1 *array, std::size_t beg, std::size_t end, T2 *second_array)
{
	std::size_t i = beg + 1;
	std::size_t j = end - 1;
	T1 m = array[beg];

	for (;;) {
		while ((i < end) && (array[i] <= m))
			i++;
		while ((j > beg) && !(array[j] <= m))
			j--;

		if (i >= j)
			break;

		std::swap(array[i], array[j]);
		std::swap(second_array[i], second_array[j]);
	}

	std::swap(array[beg], array[j]);
	std::swap(second_array[beg], second_array[j]);

	return j;
}

/**
 * version of quickSort that permutes the entries of a second array
 * according to the permutation of the first array
 * @param array array to sort
 * @param beg beginning index in array for sorting
 * @param end end-1 is the last index in array for sorting
 * @param second_array the second array is permuted according to the sort process of array
 */
template <typename T1, typename T2>
void quicksort(T1* array, std::size_t beg, std::size_t end, T2* second_array)
{
	if (beg < end) {
		std::size_t p = partition_(array, beg, end, second_array);
		quicksort(array, beg, p, second_array);
		quicksort(array, p + 1, end, second_array);
	}
}

} // end namespace BaseLib

// STL
#include <vector>

namespace BaseLib
{

template <typename T1, typename T2 = std::size_t>
void quicksortD(std::vector<T1>& array, std::size_t beg, std::size_t end, std::vector<T2>& perm)
{
	// Zip input arrays.
	std::vector<std::pair<T1, T2>> data;
	data.reserve(array.size());
	std::transform(array.begin(), array.end(), perm.begin(),
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
		array[i] = data[i].first;
		perm[i] = data[i].second;
	}
}

template <typename T1, typename T2 = std::size_t>
class Quicksort
{
public:
	Quicksort (std::vector<T1>& array, std::size_t beg, std::size_t end, std::vector<T2>& perm)
	{
		quicksort (array, beg, end, perm);
	}
private:
	std::size_t partition_(std::vector<T1>& array,
	                       std::size_t beg,
	                       std::size_t end,
	                       std::vector<T2>& perm)
	{
		std::size_t i = beg + 1;
		std::size_t j = end - 1;
		T1 m = array[beg];

		for (;;) {
			while ((i < end) && (array[i] <= m))
				i++;
			while ((j > beg) && !(array[j] <= m))
				j--;

			if (i >= j)
				break;
			std::swap(array[i], array[j]);
			std::swap(perm[i], perm[j]);
		}

		std::swap(array[beg], array[j]);
		std::swap(perm[beg], perm[j]);
		return j;
	}

	void quicksort(std::vector<T1>& array,
	               std::size_t beg,
	               std::size_t end,
	               std::vector<T2>& perm)
	{
		if (beg < end) {
			std::size_t p = partition_(array, beg, end, perm);
			quicksort(array, beg, p, perm);
			quicksort(array, p + 1, end, perm);
		}
	}
};

// specialization for pointer types
template <typename T1, typename T2>
class Quicksort <T1*, T2>
{
public:
	Quicksort (std::vector<T1*>& array, std::size_t beg, std::size_t end, std::vector<T2>& perm)
	{
		quicksort (array, beg, end, perm);
	}

	Quicksort (std::vector<std::size_t>& perm, std::size_t beg, std::size_t end,
	           std::vector<T1*>& array)
	{
		quicksort (perm, beg, end, array);
	}

private:
	std::size_t partition_(std::vector<T1*>& array,
	                       std::size_t beg,
	                       std::size_t end,
	                       std::vector<T2>& perm)
	{
		std::size_t i = beg + 1;
		std::size_t j = end - 1;
		T1* m = array[beg];

		for (;;) {
			while ((i < end) && (*array[i] <= *m))
				i++;
			while ((j > beg) && !(*array[j] <= *m))
				j--;

			if (i >= j)
				break;
			std::swap(array[i], array[j]);
			std::swap(perm[i], perm[j]);
		}

		std::swap(array[beg], array[j]);
		std::swap(perm[beg], perm[j]);
		return j;
	}

	void quicksort(std::vector<T1*>& array,
	               std::size_t beg,
	               std::size_t end,
	               std::vector<T2>& perm)
	{
		if (beg < end) {
			std::size_t p = partition_(array, beg, end, perm);
			quicksort(array, beg, p, perm);
			quicksort(array, p + 1, end, perm);
		}
	}

	std::size_t partition_(std::vector<std::size_t> &perm,
	                       std::size_t beg,
	                       std::size_t end,
	                       std::vector<T1*>& array)
	{
		std::size_t i = beg + 1;
		std::size_t j = end - 1;
		std::size_t m = perm[beg];

		for (;;) {
			while ((i < end) && (perm[i] <= m))
				i++;
			while ((j > beg) && !(perm[j] <= m))
				j--;

			if (i >= j)
				break;
			std::swap(perm[i], perm[j]);
			std::swap(array[i], array[j]);
		}

		std::swap(perm[beg], perm[j]);
		std::swap(array[beg], array[j]);
		return j;
	}

	void quicksort(std::vector<std::size_t>& perm,
	               std::size_t beg,
	               std::size_t end,
	               std::vector<T1*>& array)
	{
		if (beg < end) {
			std::size_t p = partition_(perm, beg, end, array);
			quicksort(perm, beg, p, array);
			quicksort(perm, p + 1, end, array);
		}
	}
};
} // end namespace BaseLib

#endif /* QUICKSORT_H_ */
