/**
 * \file
 * \author Thomas Fischer
 * \date   2011-02-23
 * \brief  Definition of the printList function.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PRINTLIST_H_
#define PRINTLIST_H_

// STL
#include <list>
#include <string>
#include <iostream>

namespace BaseLib {

void printList (std::list<std::size_t> const& mylist, std::string const& title)
{
	std::cout << title << std::endl;
	for (std::list<std::size_t>::const_iterator my_it (mylist.begin());
		my_it != mylist.end(); my_it++) {
		std::cout << *my_it << " ";
	}
	std::cout << std::endl;
}

} // end namespace BaseLib

#endif /* PRINTLIST_H_ */
