/*
 * printList.h
 *
 *  Created on: Feb 23, 2011
 *      Author: TF
 */

#ifndef PRINTLIST_H_
#define PRINTLIST_H_

// STL
#include <list>
#include <string>
#include <iostream>

namespace BaseLib {

void printList (std::list<size_t> const& mylist, std::string const& title)
{
	std::cout << title << std::endl;
	for (std::list<size_t>::const_iterator my_it (mylist.begin());
		my_it != mylist.end(); my_it++) {
		std::cout << *my_it << " ";
	}
	std::cout << std::endl;
}

} // end namespace BaseLib

#endif /* PRINTLIST_H_ */
