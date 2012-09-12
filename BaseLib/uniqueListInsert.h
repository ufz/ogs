/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file uniqueListInsert.h
 *
 * Created on 2011-02-23 by Thomas Fischer
 */

#ifndef UNIQUELISTINSERT_H_
#define UNIQUELISTINSERT_H_

#include <list>

namespace BaseLib {

void uniqueListInsert (std::list<std::size_t>& list, std::size_t element)
{
	// search element
	std::list<std::size_t>::const_iterator it;
	for (it = list.begin (); it != list.end(); it++) {
		if (*it == element) return;
	}
	// element not found -> insert
	list.push_back (element);
}

} // end namespace BaseLib

#endif /* UNIQUELISTINSERT_H_ */
