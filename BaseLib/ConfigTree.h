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

#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem.hpp>

extern template class boost::property_tree::basic_ptree<
    std::string, std::string, std::less<std::string>>;

namespace BaseLib
{
using ConfigTree = boost::property_tree::ptree;

boost::property_tree::ptree read_xml_config(
    boost::filesystem::path const& path);

/// Returns the JSON-representation of the given boost::property_tree.
std::string propertyTreeToString(boost::property_tree::ptree const& tree);
}

#endif  // BASELIB_CONFIGTREE_H_
