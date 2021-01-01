/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
void traverse_recursive(
    boost::property_tree::ptree& parent,
    boost::property_tree::ptree::path_type const& child_path,
    boost::property_tree::ptree& child,
    fs::path const& bench_dir,
    T& method)
{
    using boost::property_tree::ptree;

    method(parent, child_path, child, bench_dir);
    for (auto& [key, tree] : child)
    {
        ptree::path_type const cur_path = child_path / ptree::path_type(key);
        traverse_recursive(child, cur_path, tree, bench_dir, method);
    }
}

template <typename T>
void traverse(boost::property_tree::ptree& parent, const fs::path bench_dir,
              T& method)
{
    traverse_recursive(parent, "", parent, bench_dir, method);
}

void replace_includes(
    [[maybe_unused]] boost::property_tree::ptree const& parent,
    [[maybe_unused]] boost::property_tree::ptree::path_type const& child_path,
    boost::property_tree::ptree& child,
    fs::path const& bench_dir)
{
    using boost::property_tree::ptree;
    for (auto& [key, tree] : child)
    {
        if (key == "include")
        {
            auto filename = tree.get<std::string>("<xmlattr>.file");
            if (auto const filepath = fs::path(filename);
                filepath.is_relative())
            {
                filename = (bench_dir / filepath).string();
            }
            INFO("Including {:s} into project file.", filename);

            ptree include_tree;
            read_xml(filename, include_tree,
                     boost::property_tree::xml_parser::no_comments |
                         boost::property_tree::xml_parser::trim_whitespace);

            // Can only insert subtree at child
            auto& tmp_tree = child.put_child("include", include_tree);

            // Move subtree above child
            std::move(tmp_tree.begin(), tmp_tree.end(), back_inserter(child));

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
