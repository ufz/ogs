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
 * version of partition_ that additional updates the permutation vector
 * */
template <class T>
size_t partition_(T* array, size_t beg, size_t end, size_t *perm)
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

/**
 * version of quickSort that stores the permutation
 * */
template <class T>
void quicksort(T* array, size_t beg, size_t end, size_t* perm)
{
	if (beg < end) {
		size_t p = partition_(array, beg, end, perm);
		quicksort(array, beg, p, perm);
		quicksort(array, p+1, end, perm);
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
