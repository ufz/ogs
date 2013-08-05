/**
 * \file
 * \author Thomas Fischer
 * \date   2010-05-26
 * \brief  Definition of the quicksort function.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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
template <class T>
unsigned partition_(T* array, unsigned beg, unsigned end)
{
	unsigned i = beg + 1;
	unsigned j = end - 1;
	T m = array[beg];

  for (;;) {
    while ((i<end) && (array[i] < m)) i++;
    while ((j>beg) && !(array[j] < m)) j--;

    if (i >= j) break;
    std::swap(array[i], array[j]);
  }

	std::swap(array[beg], array[j]);
	return j;
}

template <class T>
void quickSort(T* array, unsigned beg, unsigned end)
{
  if (beg < end) {
    unsigned p = partition_(array, beg, end);
    quickSort(array, beg, p);
    quickSort(array, p+1, end);
  }
}

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
