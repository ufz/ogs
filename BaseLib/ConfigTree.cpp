/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConfigTree.h"

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <logog/include/logog.hpp>

// Explicitly instantiate the boost::property_tree::ptree which is a typedef to
// the following basic_ptree.
template class boost::property_tree::basic_ptree<std::string, std::string,
                                                 std::less<std::string>>;

namespace BaseLib
{

ConfigTree read_xml_config(std::string const& path)
{
	ConfigTree ptree;
	read_xml(path, ptree,
	         boost::property_tree::xml_parser::no_comments |
	             boost::property_tree::xml_parser::trim_whitespace);

	DBUG("Project configuration from file \'%s\' read.",
	     path.c_str());

	return ptree;
}

std::string propertyTreeToString(ConfigTree const& tree)
{
	std::ostringstream ss;
	boost::property_tree::write_json(ss, tree);
	return ss.str();
}

}
