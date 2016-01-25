/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef BASELIB_CONFIGTREE_H_
#define BASELIB_CONFIGTREE_H_

#include <boost/property_tree/ptree.hpp>

extern template class boost::property_tree::basic_ptree<
    std::string, std::string, std::less<std::string>>;

#endif  // BASELIB_CONFIGTREE_H_
