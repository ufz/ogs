/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConfigTreeUtil.h"

#include <boost/property_tree/xml_parser.hpp>
#include <logog/include/logog.hpp>

#include "Error.h"

namespace BaseLib
{

ConfigTreeTopLevel::ConfigTreeTopLevel(
        const std::string& filepath,
        const bool be_ruthless,
        ConfigTree::PTree&& ptree)
    : _ptree(std::move(ptree))
    , _ctree(_ptree, filepath,
             ConfigTree::onerror,
             be_ruthless ? ConfigTree::onerror : ConfigTree::onwarning)
{
}

ConfigTree const&
ConfigTreeTopLevel::operator*() const
{
    return _ctree;
}

ConfigTree const*
ConfigTreeTopLevel::operator->() const
{
    return &_ctree;
}

void
ConfigTreeTopLevel::checkAndInvalidate()
{
    ::BaseLib::checkAndInvalidate(_ctree);
}

ConfigTreeTopLevel
makeConfigTree(const std::string& filepath, const bool be_ruthless,
               const std::string& toplevel_tag)
{
    ConfigTree::PTree ptree;

    // note: Trimming whitespace and ignoring comments is crucial in order
    //       for our configuration tree implementation to work!
    try
    {
        read_xml(filepath, ptree,
                 boost::property_tree::xml_parser::no_comments |
                     boost::property_tree::xml_parser::trim_whitespace);
    }
    catch (boost::property_tree::xml_parser_error const& e)
    {
        OGS_FATAL("Error while parsing XML file `%s' at line %lu: %s.",
                  e.filename().c_str(), e.line(), e.message().c_str());
    }

    DBUG("Project configuration from file \'%s\' read.", filepath.c_str());

    if (auto child = ptree.get_child_optional(toplevel_tag)) {
        return ConfigTreeTopLevel(filepath, be_ruthless, std::move(*child));
    }
    OGS_FATAL("Tag <%s> has not been found in file `%s'.", toplevel_tag.c_str(),
              filepath.c_str());
}

}
