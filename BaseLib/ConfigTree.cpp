/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConfigTree.h"

// Explicitly instantiate the boost::property_tree::ptree which is a typedef to
// the following basic_ptree.
template class boost::property_tree::basic_ptree<std::string, std::string,
                                                 std::less<std::string>>;
