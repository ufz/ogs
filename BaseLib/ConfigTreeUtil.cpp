/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConfigTreeUtil.h"

#include <boost/property_tree/xml_parser.hpp>
#include <regex>

#include "Error.h"
#include "Logging.h"

namespace BaseLib
{
ConfigTree makeConfigTree(const std::string& filepath, const bool be_ruthless,
                          const std::string& toplevel_tag,
                          std::stringstream& prj_stream)
{
    ConfigTree::PTree ptree;

    // note: Trimming whitespace and ignoring comments is crucial in order
    //       for our configuration tree implementation to work!
    try
    {
        read_xml(prj_stream, ptree,
                 boost::property_tree::xml_parser::no_comments |
                     boost::property_tree::xml_parser::trim_whitespace);
    }
    catch (boost::property_tree::xml_parser_error const& e)
    {
        OGS_FATAL("Error while parsing XML file `{:s}' at line {:d}: {:s}.",
                  e.filename(), e.line(), e.message());
    }

    DBUG("Project configuration from file '{:s}' read.", filepath);

    if (auto opt_child = ptree.get_child_optional(toplevel_tag))
    {
        auto const callback =
            be_ruthless ? ConfigTree::onerror : ConfigTree::onwarning;

        return ConfigTree(std::move(*opt_child), filepath, callback, callback);
    }
    OGS_FATAL("Tag <{:s}> has not been found in file `{:s}'.", toplevel_tag,
              filepath);
}

ConfigTree makeConfigTreeFromFile(const std::string& filepath,
                                  const bool be_ruthless,
                                  const std::string& toplevel_tag)
{
    std::ifstream file(filepath);
    if (file)
    {
        std::stringstream buffer;
        buffer << file.rdbuf();
        file.close();

        return makeConfigTree(filepath, be_ruthless, toplevel_tag, buffer);
    }
    else
    {
        OGS_FATAL("Could not read from file {:s}!", filepath);
    }
}

}  // namespace BaseLib
