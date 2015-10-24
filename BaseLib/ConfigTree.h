/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef BASELIB_CONFIGTREE_H_
#define BASELIB_CONFIGTREE_H_

#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/filesystem.hpp>

namespace BaseLib
{

boost::property_tree::ptree read_xml_config(
    boost::filesystem::path const& path);
}

#endif  // BASELIB_CONFIGTREE_H_
