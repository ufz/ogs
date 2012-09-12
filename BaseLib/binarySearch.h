/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file binarySearch.h
 *
 * Created on 2010-06-07 by Thomas Fischer
 *
 */

// STL
#include <limits>
#include <vector>
#include <cstddef>

#ifndef BINARYSEARCH_H_
#define BINARYSEARCH_H_

namespace BaseLib {

/**
 * Binary search in a sorted vector of elements to get the
 * id of an element according its key.
 * @param key the key for the element
 * @param beg beginning index in the sorted vector of elements
 * @param end ending index in the sorted vector of elements
 * @param array the vector of elements
 * @return the id of the element in the vector or, if not found,
 * the value std::numeric_limits<std::size_t>::max()
 */
template <class T>
std::size_t searchElement (const T& key, std::size_t beg, std::size_t end, const std::vector<T>& array)
{
	if (beg >= end) return std::numeric_limits<std::size_t>::max();
	std::size_t m ((end+beg)/2);

	if (key == array[m]) {
		return m;
	}
	if (key < array[m]) {
		return searchElement (key, beg, m, array);
	}
	return searchElement (key, m+1, end, array);
}

std::size_t searchElement (double const& val, std::size_t beg, std::size_t end, const std::vector<double>& array);

} // end namespace BaseLib

#endif /* BINARYSEARCH_H_ */
