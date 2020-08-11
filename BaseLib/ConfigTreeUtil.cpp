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

#include <boost/bind.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <regex>

#include "Error.h"
#include "Logging.h"
#include "filesystem.h"

namespace BaseLib
{
ConfigTreeTopLevel::ConfigTreeTopLevel(const std::string& filepath,
                                       const bool be_ruthless,
                                       ConfigTree::PTree&& ptree)
    : _ptree(std::move(ptree)),
      _ctree(_ptree, filepath, ConfigTree::onerror,
             be_ruthless ? ConfigTree::onerror : ConfigTree::onwarning)
{
}

ConfigTree const& ConfigTreeTopLevel::operator*() const
{
    return _ctree;
}

ConfigTree const* ConfigTreeTopLevel::operator->() const
{
    return &_ctree;
}

void ConfigTreeTopLevel::checkAndInvalidate()
{
    ::BaseLib::checkAndInvalidate(_ctree);
}

// Adapted from
// https://stackoverflow.com/questions/8154107/how-do-i-merge-update-a-boostproperty-treeptree/8175833
template <typename T>
void traverse_recursive(boost::property_tree::ptree& parent,
                        const boost::property_tree::ptree::path_type& childPath,
                        boost::property_tree::ptree& child,
                        const fs::path benchDir,
                        T& method)
{
    using boost::property_tree::ptree;

    method(parent, childPath, child, benchDir);
    for (ptree::iterator it = child.begin(); it != child.end(); ++it)
    {
        ptree::path_type curPath = childPath / ptree::path_type(it->first);
        traverse_recursive(child, curPath, it->second, benchDir, method);
    }
}

template <typename T>
void traverse(boost::property_tree::ptree& parent, const fs::path benchDir,
              T& method)
{
    traverse_recursive(parent, "", parent, benchDir, method);
}

void replace_includes(
    [[maybe_unused]] const boost::property_tree::ptree& parent,
    [[maybe_unused]] const boost::property_tree::ptree::path_type& childPath,
    boost::property_tree::ptree& child,
    const fs::path benchDir)
{
    using boost::property_tree::ptree;
    for (ptree::const_iterator it = child.begin(); it != child.end(); ++it)
    {
        if (it->first == "include")
        {
            auto filename = it->second.get<std::string>("<xmlattr>.file");
            auto filepath = fs::path(filename);
            if (filepath.is_relative())
            {
                filename = (benchDir / filepath).string();
            }
            INFO("Including {:s} into project file.", filename);

            ptree includeTree;
            read_xml(filename, includeTree,
                     boost::property_tree::xml_parser::no_comments |
                         boost::property_tree::xml_parser::trim_whitespace);

            // Can only insert subtree at child
            auto& tmpTree = child.put_child("include", includeTree);

            // Move subtree above child
            std::move(tmpTree.begin(), tmpTree.end(), back_inserter(child));

            // Erase child
            child.erase("include");

            // There can only be one include under a parent element!
            break;
        }
    }
}

ConfigTreeTopLevel makeConfigTree(const std::string& filepath,
                                  const bool be_ruthless,
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

        if (toplevel_tag == "OpenGeoSysProject")
        {
            traverse(ptree, fs::path(filepath).parent_path(), replace_includes);
        }
    }
    catch (boost::property_tree::xml_parser_error const& e)
    {
        OGS_FATAL("Error while parsing XML file `{:s}' at line {:d}: {:s}.",
                  e.filename(), e.line(), e.message());
    }

    DBUG("Project configuration from file '{:s}' read.", filepath);

    if (auto child = ptree.get_child_optional(toplevel_tag))
    {
        return ConfigTreeTopLevel(filepath, be_ruthless, std::move(*child));
    }
    OGS_FATAL("Tag <{:s}> has not been found in file `{:s}'.", toplevel_tag,
              filepath);
}

}  // namespace BaseLib
