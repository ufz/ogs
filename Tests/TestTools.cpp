// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Tests/TestTools.h"

#include <boost/property_tree/xml_parser.hpp>

namespace Tests
{
boost::property_tree::ptree readXml(const char xml[])
{
    boost::property_tree::ptree ptree;
    std::istringstream xml_str(xml);
    boost::property_tree::read_xml(
        xml_str, ptree,
        boost::property_tree::xml_parser::no_comments |
            boost::property_tree::xml_parser::trim_whitespace);
    return ptree;
}

}  // namespace Tests
