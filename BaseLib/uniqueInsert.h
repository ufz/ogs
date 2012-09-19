/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file uniqueListInsert.h
 *
 * Created on 2011-02-23 by Thomas Fischer
 */

#ifndef UNIQUELISTINSERT_H_
#define UNIQUELISTINSERT_H_

#include <list>
#include <vector>

namespace BaseLib {

template<typename T>
void uniqueListInsert (std::list<T>& list, T element)
{
	// search element
	typename std::list<T>::const_iterator it;
	for (it = list.begin (); it != list.end(); it++) {
		if (*it == element) return;
	}
	// element not found -> insert
	list.push_back (element);
}

template<typename T>
void uniqueVectorInsert (std::vector<T>& vec, T element)
{
	// search element
	typename std::vector<T>::const_iterator it;
	for (it = vec.begin (); it != vec.end(); ++it)
		if (*it == element) return;
	// element not found -> insert
	vec.push_back (element);
}


} // end namespace BaseLib

#endif /* UNIQUELISTINSERT_H_ */
