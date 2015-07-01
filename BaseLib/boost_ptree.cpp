/*
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include <boost/property_tree/ptree.hpp>

template class boost::property_tree::basic_ptree<std::basic_string<char>,
      std::basic_string<char>, std::less<std::basic_string<char> > >;
