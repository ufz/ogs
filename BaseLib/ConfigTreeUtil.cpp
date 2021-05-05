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

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <spdlog/spdlog.h>
#include <xml_patch.h>

#include <boost/property_tree/xml_parser.hpp>
#include <regex>

#include "BaseLib/FileTools.h"
#include "Error.h"
#include "Logging.h"
#include "filesystem.h"

namespace BaseLib
{
ConfigTreeTopLevel::ConfigTreeTopLevel(const std::string& filepath,
                                       const bool be_ruthless,
                                       ConfigTree::PTree&& ptree)
    : ptree_(std::move(ptree)),
      ctree_(ptree_, filepath, ConfigTree::onerror,
             be_ruthless ? ConfigTree::onerror : ConfigTree::onwarning)
{
}

ConfigTree const& ConfigTreeTopLevel::operator*() const
{
    return ctree_;
}

ConfigTree const* ConfigTreeTopLevel::operator->() const
{
    return &ctree_;
}

void ConfigTreeTopLevel::checkAndInvalidate()
{
    ::BaseLib::checkAndInvalidate(ctree_);
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
                                  const std::string& toplevel_tag,
                                  std::vector<std::string>
                                      patch_files,
                                  std::string outdir)
{
    std::string prj_file = filepath;
    std::string includepath = "";
    if (BaseLib::getFileExtension(prj_file) == ".xml")
    {
        if (!patch_files.empty())
        {
            OGS_FATAL(
                "It is not allowed to specify additional patch files "
                "if a patch file was already specified as the "
                "prj-file.");
        }
        auto patch = xmlParseFile(prj_file.c_str());
        auto node = xmlDocGetRootElement(patch);
        auto base_file = xmlGetProp(node, (const xmlChar*)"base_file");
        if (base_file == nullptr)
        {
            OGS_FATAL("Error reading base prj file in given patch file {:s}.",
                      prj_file);
        }
        patch_files = {prj_file};
        std::stringstream ss;
        ss << base_file;
        prj_file = BaseLib::joinPaths(BaseLib::extractPath(prj_file), ss.str());
    }

    if (!patch_files.empty())
    {
        includepath = fs::canonical(fs::path(prj_file)).parent_path().string();
        std::string current_prj_file = prj_file;
        std::string current_prj_file_base =
            BaseLib::extractBaseNameWithoutExtension(current_prj_file);
        for (const auto& patch_file : patch_files)
        {
            auto patch = xmlParseFile(patch_file.c_str());
            if (patch == NULL)
            {
                OGS_FATAL("Error reading XML diff file {:s}.", patch_file);
            }

            auto doc = xmlParseFile(current_prj_file.c_str());
            if (doc == NULL)
            {
                OGS_FATAL("Error reading project file {:s}.", current_prj_file);
            }

            auto node = xmlDocGetRootElement(patch);
            int rc = 0;
            for (node = node ? node->children : NULL; node; node = node->next)
            {
                if (node->type != XML_ELEMENT_NODE)
                    continue;

                if (!strcmp((char*)node->name, "add"))
                    rc = xml_patch_add(doc, node);
                else if (!strcmp((char*)node->name, "replace"))
                    rc = xml_patch_replace(doc, node);
                else if (!strcmp((char*)node->name, "remove"))
                    rc = xml_patch_remove(doc, node);
                else
                    continue;

                if (rc)
                {
                    OGS_FATAL(
                        "Error while patching prj file {:s} with patch "
                        "file "
                        "{:}.",
                        filepath, patch_file);
                }
            }

            current_prj_file_base =
                current_prj_file_base + "_" +
                BaseLib::extractBaseNameWithoutExtension(patch_file);
            current_prj_file =
                BaseLib::joinPaths(outdir, current_prj_file_base + ".prj");

            // TODO: make name unique
            xmlSaveFile(current_prj_file.c_str(), doc);

            xmlFreeDoc(doc);
            xmlFreeDoc(patch);
        }
        xmlCleanupParser();
        prj_file = current_prj_file;
    }
    else
    {
        prj_file = filepath;
    }

    ConfigTree::PTree ptree;

    // note: Trimming whitespace and ignoring comments is crucial in order
    //       for our configuration tree implementation to work!
    try
    {
        read_xml(prj_file, ptree,
                 boost::property_tree::xml_parser::no_comments |
                     boost::property_tree::xml_parser::trim_whitespace);

        if (toplevel_tag == "OpenGeoSysProject")
        {
            if (includepath == "")
            {
                includepath = fs::path(filepath).parent_path().string();
            }
            traverse(ptree, fs::path(includepath), replace_includes);
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
