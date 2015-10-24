/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConfigTree.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <logog/include/logog.hpp>

namespace BaseLib
{
boost::property_tree::ptree read_xml_config(boost::filesystem::path const& path)
{
	boost::property_tree::ptree ptree;
	read_xml(path.string(), ptree,
	         boost::property_tree::xml_parser::no_comments |
	             boost::property_tree::xml_parser::trim_whitespace);

	DBUG("Project configuration from file \'%s\' read.",
	     path.string().c_str());

	return ptree;
}
}
