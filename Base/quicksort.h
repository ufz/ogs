/*
 * quicksort.h
 *
 *  Created on: May 26, 2010
 *      Author: TF
 */

#ifndef QUICKSORT_H_
#define QUICKSORT_H_

// STL
#include <cstddef>

// Base
#include "swap.h"

template <class T>
void quickSort(T* array, unsigned beg, unsigned end)
{
  if (beg < end) {
    unsigned p = partition_(array, beg, end);
    quickSort(array, beg, p);
    quickSort(array, p+1, end);
  }
}

template <class T>
unsigned partition_(T* array, unsigned beg, unsigned end)
{
  unsigned i = beg+1;
  unsigned j = end-1;
  T m = array[beg];

  for (;;) {
    while ((i<end) && (array[i] < m)) i++;
    while ((j>beg) && !(array[j] < m)) j--;

    if (i >= j) break;
    BaseLib::swap(array[i], array[j]);
  }

  BaseLib::swap(array[beg], array[j]);
  return j;
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
 * @return
 */
template <typename T1, typename T2>
size_t partition_(T1* array, size_t beg, size_t end, T2 *second_array)
{
	size_t i = beg + 1;
	size_t j = end - 1;
	T1 m = array[beg];

	for (;;) {
		while ((i < end) && (array[i] <= m))
			i++;
		while ((j > beg) && !(array[j] <= m))
			j--;

		if (i >= j) break;

		BaseLib::swap(array[i], array[j]);
		BaseLib::swap(second_array[i], second_array[j]);
	}

	BaseLib::swap(array[beg], array[j]);
	BaseLib::swap(second_array[beg], second_array[j]);

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
void quicksort(T1* array, size_t beg, size_t end, T2* second_array)
{
	if (beg < end) {
		size_t p = partition_(array, beg, end, second_array);
		quicksort(array, beg, p, second_array);
		quicksort(array, p+1, end, second_array);
	}
}

// STL
#include <vector>

template <typename T>
class Quicksort {
public:
	Quicksort (std::vector<T>& array, size_t beg, size_t end, std::vector<size_t>& perm)
	{
		quicksort (array, beg, end, perm);
	}
private:
	size_t partition_(std::vector<T>& array, size_t beg, size_t end, std::vector<size_t>& perm)
	{
		size_t i = beg + 1;
		size_t j = end - 1;
		T m = array[beg];

		for (;;) {
			while ((i < end) && (array[i] <= m))
				i++;
			while ((j > beg) && !(array[j] <= m))
				j--;

			if (i >= j)
				break;
			BaseLib::swap(array[i], array[j]);
			BaseLib::swap(perm[i], perm[j]);
		}

		BaseLib::swap(array[beg], array[j]);
		BaseLib::swap(perm[beg], perm[j]);
		return j;
	}

	void quicksort(std::vector<T>& array, size_t beg, size_t end, std::vector<size_t>& perm)
	{
		if (beg < end) {
			size_t p = partition_(array, beg, end, perm);
			quicksort(array, beg, p, perm);
			quicksort(array, p+1, end, perm);
		}
	}
};

// specialization for pointer types
template <typename T>
class Quicksort <T *> {
public:
	Quicksort (std::vector<T*>& array, size_t beg, size_t end, std::vector<size_t>& perm)
	{
		quicksort (array, beg, end, perm);
	}

	Quicksort (std::vector<size_t>& perm, size_t beg, size_t end, std::vector<T*>& array)
	{
		quicksort (perm, beg, end, array);
	}

private:
	size_t partition_(std::vector<T*>& array, size_t beg, size_t end, std::vector<size_t>& perm)
	{
		size_t i = beg + 1;
		size_t j = end - 1;
		T* m = array[beg];

		for (;;) {
			while ((i < end) && (*array[i] <= *m))
				i++;
			while ((j > beg) && !(*array[j] <= *m))
				j--;

			if (i >= j)
				break;
			BaseLib::swap(array[i], array[j]);
			BaseLib::swap(perm[i], perm[j]);
		}

		BaseLib::swap(array[beg], array[j]);
		BaseLib::swap(perm[beg], perm[j]);
		return j;
	}

	void quicksort(std::vector<T*>& array, size_t beg, size_t end, std::vector<size_t>& perm)
	{
		if (beg < end) {
			size_t p = partition_(array, beg, end, perm);
			quicksort(array, beg, p, perm);
			quicksort(array, p+1, end, perm);
		}
	}

	size_t partition_(std::vector<size_t> &perm, size_t beg, size_t end, std::vector<T*>& array)
	{
		size_t i = beg + 1;
		size_t j = end - 1;
		size_t m = perm[beg];

		for (;;) {
			while ((i < end) && (perm[i] <= m))
				i++;
			while ((j > beg) && !(perm[j] <= m))
				j--;

			if (i >= j)
				break;
			BaseLib::swap(perm[i], perm[j]);
			BaseLib::swap(array[i], array[j]);
		}

		BaseLib::swap(perm[beg], perm[j]);
		BaseLib::swap(array[beg], array[j]);
		return j;
	}

	void quicksort(std::vector<size_t>& perm, size_t beg, size_t end, std::vector<T*>& array)
	{
		if (beg < end) {
			size_t p = partition_(perm, beg, end, array);
			quicksort(perm, beg, p, array);
			quicksort(perm, p+1, end, array);
		}
	}

};

#endif /* QUICKSORT_H_ */
