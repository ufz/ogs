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

#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem.hpp>

namespace BaseLib
{
using ConfigTree = boost::property_tree::ptree;

boost::property_tree::ptree read_xml_config(
    boost::filesystem::path const& path);

/// Returns the JSON-representation of the given boost::property_tree.
std::string propertyTreeToString(boost::property_tree::ptree const& tree);
}

#endif  // BASELIB_CONFIGTREE_H_
