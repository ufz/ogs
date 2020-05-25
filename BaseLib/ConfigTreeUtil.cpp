/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConfigTreeUtil.h"

#include <boost/property_tree/xml_parser.hpp>
#include "Logging.h"

#include "Error.h"

namespace BaseLib
{

ConfigTreeTopLevel::ConfigTreeTopLevel(
        const std::string& filepath,
        const bool be_ruthless,
        ConfigTree::PTree&& ptree)
    : ptree_(std::move(ptree))
    , ctree_(ptree_, filepath,
             ConfigTree::onerror,
             be_ruthless ? ConfigTree::onerror : ConfigTree::onwarning)
{
}

ConfigTree const&
ConfigTreeTopLevel::operator*() const
{
    return ctree_;
}

ConfigTree const*
ConfigTreeTopLevel::operator->() const
{
    return &ctree_;
}

void
ConfigTreeTopLevel::checkAndInvalidate()
{
    ::BaseLib::checkAndInvalidate(ctree_);
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
        OGS_FATAL("Error while parsing XML file `{:s}' at line {:d}: {:s}.",
                  e.filename(), e.line(), e.message());
    }

    DBUG("Project configuration from file '{:s}' read.", filepath);

    if (auto child = ptree.get_child_optional(toplevel_tag)) {
        return ConfigTreeTopLevel(filepath, be_ruthless, std::move(*child));
    }
    OGS_FATAL("Tag <{:s}> has not been found in file `{:s}'.", toplevel_tag,
              filepath);
}

}  // namespace BaseLib
