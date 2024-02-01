/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PrjProcessing.h"

#include <libxml/parser.h>
#include <libxml/xmlstring.h>
#include <spdlog/fmt/bundled/core.h>
#include <xml_patch.h>

#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>

#include "DisableFPE.h"
#include "Error.h"
#include "FileTools.h"
#include "Logging.h"

namespace
{
std::string iostateToString(std::ios_base::iostate const state)
{
    std::string result;

    if (state == std::ios_base::goodbit)
    {
        result = "goodbit";
    }
    else
    {
        if (state & std::ios_base::eofbit)
        {
            result += "eofbit ";
        }
        if (state & std::ios_base::failbit)
        {
            result += "failbit ";
        }
        if (state & std::ios_base::badbit)
        {
            result += "badbit";
        }
        // Remove trailing space if there is one
        if (!result.empty() && result.back() == ' ')
        {
            result.pop_back();
        }
    }
    return result;
}
}  // namespace

namespace BaseLib
{
void traverseIncludes(xmlDoc* doc, xmlNode* node,
                      std::filesystem::path const& prj_dir)
{
    xmlNode* cur_node = nullptr;
    for (cur_node = node; cur_node; cur_node = cur_node->next)
    {
        if (cur_node->type != XML_ELEMENT_NODE)
        {
            continue;
        }
        if (xmlStrEqual(cur_node->name, xmlCharStrdup("include")))
        {
            auto include_file_char_pointer =
                xmlGetProp(cur_node, xmlCharStrdup("file"));
            if (include_file_char_pointer == nullptr)
            {
                OGS_FATAL(
                    "Error while processing includes in prj file. Error in "
                    "element '{:s}' on line {:d}: no file attribute given!",
                    reinterpret_cast<const char*>(cur_node->name),
                    cur_node->line);
            }
            auto filename_length = xmlStrlen(include_file_char_pointer);
            std::string filename(
                reinterpret_cast<char*>(include_file_char_pointer),
                filename_length);
            if (auto const filepath = std::filesystem::path(filename);
                filepath.is_relative())
            {
                filename = (prj_dir / filepath).string();
            }

            if (!std::filesystem::exists(filename))
            {
                OGS_FATAL(
                    "Error while processing includes in prj file. Error in "
                    "element '{:s}' on line {:d}: Include file is not "
                    "existing: "
                    "{:s}!",
                    reinterpret_cast<const char*>(cur_node->name),
                    cur_node->line,
                    reinterpret_cast<const char*>(include_file_char_pointer));
            }
            INFO("Including {:s} into project file.", filename);

            const std::ifstream input_stream(filename, std::ios_base::binary);
            if (input_stream.fail())
            {
                OGS_FATAL("Failed to open file {}!", filename);
            }
            std::stringstream buffer;
            buffer << input_stream.rdbuf();
            const std::string xml = buffer.str();

            // Replace lines containing <?xml ... ?>
            std::regex xmlDeclaration("<\\?xml.*?>");
            const std::string xml_filtered =
                std::regex_replace(xml, xmlDeclaration, "");

            xmlNodePtr pNewNode = nullptr;
            xmlParseInNodeContext(cur_node->parent, xml_filtered.c_str(),
                                  (int)xml_filtered.length(), 0, &pNewNode);
            if (pNewNode != nullptr)
            {
                // add new xml node to parent
                xmlNode* pChild = pNewNode;
                while (pChild != nullptr)
                {
                    xmlAddChild(cur_node->parent, xmlCopyNode(pChild, 1));
                    pChild = pChild->next;
                }
                xmlFreeNode(pNewNode);
            }

            // Cleanup and continue on next node
            auto next_node = cur_node->next;
            xmlUnlinkNode(cur_node);
            xmlFreeNode(cur_node);
            cur_node = next_node;
        }
        traverseIncludes(doc, cur_node->children, prj_dir);
    }
    xmlFreeNode(cur_node);
}

void replaceIncludes(std::stringstream& prj_stream,
                     std::filesystem::path const& prj_dir)
{
    // Parsing the XML triggers floating point exceptions. Because we are not
    // debugging libxml2 (or other libraries) at this point, the floating point
    // exceptions are temporarily disabled and are restored at the end of the
    // function.
    [[maybe_unused]] DisableFPE disable_fpe;

    auto doc = xmlReadMemory(prj_stream.str().c_str(), prj_stream.str().size(),
                             nullptr, nullptr, 0);
    if (doc == nullptr)
    {
        OGS_FATAL("Error reading project file from memory.");
    }

    auto root_node = xmlDocGetRootElement(doc);
    traverseIncludes(doc, root_node, prj_dir);

    xmlChar* xmlbuff;
    int buffersize;
    xmlDocDumpMemory(doc, &xmlbuff, &buffersize);
    prj_stream.str("");  // Empty stream
    prj_stream << xmlbuff;

    xmlFree(xmlbuff);
    xmlFreeDoc(doc);
}

// Applies a patch file to the prj content in prj_stream.
void patchStream(std::string const& patch_file, std::stringstream& prj_stream,
                 bool after_includes = false)
{
    // xmlReadFile(patch_file.c_str(), nullptr, 0); leads to malloc errors in
    // xmlFreeDoc(doc) below
    auto patch = xmlParseFile(patch_file.c_str());
    if (patch == nullptr)
    {
        OGS_FATAL("Error reading XML diff file {:s}.", patch_file);
    }

    auto doc = xmlReadMemory(prj_stream.str().c_str(), prj_stream.str().size(),
                             nullptr, nullptr, 0);
    if (doc == nullptr)
    {
        OGS_FATAL("Error reading project file from memory.");
    }

    auto node = xmlDocGetRootElement(patch);
    int rc = 0;
    for (node = node ? node->children : nullptr; node; node = node->next)
    {
        if (node->type != XML_ELEMENT_NODE)
        {
            continue;
        }
        bool node_after_includes = false;
        xmlAttr* attribute = node->properties;
        while (attribute)
        {
            // Check for after_includes-attribute
            xmlChar* value =
                xmlNodeListGetString(node->doc, attribute->children, 1);
            if (xmlStrEqual(attribute->name, xmlCharStrdup("after_includes")) &&
                xmlStrEqual(value, xmlCharStrdup("true")))
            {
                node_after_includes = true;
            }

            xmlFree(value);
            attribute = attribute->next;
        }

        if (after_includes != node_after_includes)
        {
            continue;
        }

        if (xmlStrEqual(node->name, xmlCharStrdup("add")))
        {
            rc = xml_patch_add(doc, node);
        }
        else if (xmlStrEqual(node->name, xmlCharStrdup("replace")))
        {
            rc = xml_patch_replace(doc, node);
        }
        else if (xmlStrEqual(node->name, xmlCharStrdup("remove")))
        {
            rc = xml_patch_remove(doc, node);
        }
        else
        {
            OGS_FATAL(
                "Error while patching prj file with patch file {:}. Only "
                "'add', 'replace' and 'remove' elements are allowed! Got an "
                "element '{:s}' on line {:d}.",
                patch_file, reinterpret_cast<const char*>(node->name),
                node->line);
        }

        if (rc)
        {
            OGS_FATAL(
                "Error while patching prj file with patch file {:}. Error in "
                "element '{:s}' on line {:d}.",
                patch_file, reinterpret_cast<const char*>(node->name),
                node->line);
        }
    }

    xmlChar* xmlbuff;
    int buffersize;
    xmlDocDumpMemory(doc, &xmlbuff, &buffersize);
    prj_stream.str("");  // Empty stream
    prj_stream << xmlbuff;

    xmlFree(xmlbuff);
    xmlFreeDoc(doc);
    xmlFreeDoc(patch);
}

// Will set prj_file to the actual .prj file and returns the final prj file
// content in prj_stream.
void readAndPatchPrj(std::stringstream& prj_stream, std::string& prj_file,
                     std::vector<std::string>& patch_files)
{
    // Extract base project file path if an xml (patch) file is given as prj
    // file and it contains the base_file attribute.
    if (BaseLib::getFileExtension(prj_file) == ".xml")
    {
        if (!patch_files.empty())
        {
            OGS_FATAL(
                "It is not allowed to specify additional patch files "
                "if a patch file was already specified as the "
                "prj-file.");
        }
        auto patch = xmlReadFile(prj_file.c_str(), nullptr, 0);
        auto node = xmlDocGetRootElement(patch);
        xmlChar const base_file_string[] = "base_file";
        auto base_file = xmlGetProp(node, base_file_string);
        if (base_file == nullptr)
        {
            OGS_FATAL(
                "Error reading base prj file (base_file attribute) in given "
                "patch file {:s}.",
                prj_file);
        }
        patch_files = {prj_file};
        std::stringstream ss;
        ss << base_file;
        prj_file = BaseLib::joinPaths(BaseLib::extractPath(prj_file), ss.str());
    }

    // read base prj file into stream
    if (std::ifstream file(prj_file); file)
    {
        prj_stream << file.rdbuf();
    }
    else
    {
        if (!BaseLib::IsFileExisting(prj_file))
        {
            ERR("File {:s} does not exist.", prj_file);
        }
        DBUG("Stream state flags: {:s}.", iostateToString(file.rdstate()));
        OGS_FATAL("Could not open project file '{:s}' for reading.", prj_file);
    }

    // apply xml patches to stream
    for (const auto& patch_file : patch_files)
    {
        patchStream(patch_file, prj_stream);
    }
}

void prepareProjectFile(std::stringstream& prj_stream,
                        const std::string& filepath,
                        const std::vector<std::string>& patch_files,
                        bool write_prj, const std::string& out_directory)
{
    std::string prj_file = filepath;

    std::vector<std::string> patch_files_copy = patch_files;
    readAndPatchPrj(prj_stream, prj_file, patch_files_copy);
    replaceIncludes(prj_stream,
                    std::filesystem::absolute(std::filesystem::path(prj_file))
                        .parent_path());
    // re-apply xml patches to stream
    for (const auto& patch_file : patch_files_copy)
    {
        patchStream(patch_file, prj_stream, true);
    }

    if (write_prj)
    {
        // The following two lines should set indenting to 4 spaces but it does
        // not work. 2 spaces are default.
        //
        // xmlThrDefIndentTreeOutput(1);
        // xmlThrDefTreeIndentString("    "); // 4 spaces indent
        // XML_PARSE_NOBLANKS -> pretty-print
        auto doc =
            xmlReadMemory(prj_stream.str().c_str(), prj_stream.str().size(),
                          nullptr, nullptr, XML_PARSE_NOBLANKS);
        auto prj_out = (std::filesystem::path(out_directory) /
                        std::filesystem::path(filepath).stem())
                           .string() +
                       "_processed.prj";
        xmlSaveFormatFileEnc(prj_out.c_str(), doc, "utf-8", 1);
        INFO("Processed project file written to {:s}.", prj_out);
        xmlFreeDoc(doc);
    }
    xmlCleanupParser();
}
}  // namespace BaseLib
