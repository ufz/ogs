/**
 * \file uniqueListInsert.h
 *
 * Created on 2011-02-23 by Thomas Fischer
 */

#ifndef UNIQUELISTINSERT_H_
#define UNIQUELISTINSERT_H_

#include <list>

namespace BaseLib {

void uniqueListInsert (std::list<size_t>& list, size_t element)
{
	// search element
	std::list<size_t>::const_iterator it;
	for (it = list.begin (); it != list.end(); it++) {
		if (*it == element) return;
	}
	// element not found -> insert
	list.push_back (element);
}

} // end namespace BaseLib

#endif /* UNIQUELISTINSERT_H_ */
