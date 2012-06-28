/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com/LICENSE.txt
 *
 *
 * \file printList.h
 *
 * Created on 2011-02-23 by Thomas Fischer
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
