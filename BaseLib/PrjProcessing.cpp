/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PrjProcessing.h"

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <xml_patch.h>

#include "Error.h"
#include "Logging.h"

namespace BaseLib
{

bool isEqual(const unsigned char* a, std::string const& b)
{
    return std::string{reinterpret_cast<const char*>(a)} == b;
}

void traverseIncludes(xmlDoc* doc, xmlNode* node,
                      std::filesystem::path const& prj_dir)
{
    xmlNode* cur_node = nullptr;
    for (cur_node = node; cur_node; cur_node = cur_node->next)
    {
        if (cur_node->type != XML_ELEMENT_NODE ||
            cur_node->type == XML_TEXT_NODE)
        {
            continue;
        }
        if (isEqual(cur_node->name, "include"))
        // if (!strcmp(reinterpret_cast<const char*>(cur_node->name),
        // "include"))
        {
            auto include_file_char_pointer =
                xmlGetProp(cur_node, xmlCharStrdup("file"));
            if (include_file_char_pointer == nullptr)
            {
                OGS_FATAL(
                    "Error while processing includes in prj file. Error in "
                    "element '{:s}' on line {:d}: no file attribute given!",
                    cur_node->name, cur_node->line);
            }
            std::string filename(
                reinterpret_cast<char*>(include_file_char_pointer));
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
                    cur_node->name, cur_node->line, include_file_char_pointer);
            }
            INFO("Including {:s} into project file.", filename);

            const std::ifstream input_stream(filename, std::ios_base::binary);
            if (input_stream.fail())
            {
                OGS_FATAL("Failed to open file {s}!", filename);
            }
            std::stringstream buffer;
            buffer << input_stream.rdbuf();
            const std::string xml = buffer.str();

            xmlNodePtr pNewNode = nullptr;
            xmlParseInNodeContext(cur_node->parent, xml.c_str(),
                                  (int)xml.length(), 0, &pNewNode);
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
    auto doc =
        xmlParseMemory(prj_stream.str().c_str(), prj_stream.str().size());
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
    auto patch = xmlParseFile(patch_file.c_str());
    if (patch == nullptr)
    {
        OGS_FATAL("Error reading XML diff file {:s}.", patch_file);
    }

    auto doc =
        xmlParseMemory(prj_stream.str().c_str(), prj_stream.str().size());
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
            if (isEqual(attribute->name, "after_includes") &&
                isEqual(value, "true"))
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

        if (isEqual(node->name, "add"))
        {
            rc = xml_patch_add(doc, node);
        }
        else if (isEqual(node->name, "replace"))
        {
            rc = xml_patch_replace(doc, node);
        }
        else if (isEqual(node->name, "remove"))
        {
            rc = xml_patch_remove(doc, node);
        }
        else
        {
            OGS_FATAL(
                "Error while patching prj file with patch file {:}. Only "
                "'add', 'replace' and 'remove' elements are allowed! Got an "
                "element '{:s}' on line {:d}.",
                patch_file, node->name, node->line);
        }

        if (rc)
        {
            OGS_FATAL(
                "Error while patching prj file with patch file {:}. Error in "
                "element '{:s}' on line {:d}.",
                patch_file, node->name, node->line);
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
        auto patch = xmlParseFile(prj_file.c_str());
        auto node = xmlDocGetRootElement(patch);
        auto base_file = xmlGetProp(node, (const xmlChar*)"base_file");
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
        DBUG("Stream state flags: {:b}.", file.rdstate());
        OGS_FATAL("Could not open project file '{:s}' for reading.", prj_file);
    }

    // apply xml patches to stream
    if (!patch_files.empty())
    {
        for (const auto& patch_file : patch_files)
        {
            patchStream(patch_file, prj_stream);
        }
    }
}

void prepareProjectFile(std::stringstream& prj_stream,
                        const std::string& filepath,
                        const std::vector<std::string>& patch_files)
{
    std::string prj_file = filepath;

    std::vector<std::string> patch_files_copy = patch_files;
    readAndPatchPrj(prj_stream, prj_file, patch_files_copy);
    replaceIncludes(prj_stream, std::filesystem::path(prj_file).parent_path());
    // re-apply xml patches to stream
    if (!patch_files_copy.empty())
    {
        for (const auto& patch_file : patch_files_copy)
        {
            patchStream(patch_file, prj_stream, true);
        }
    }
    if (true)
    {
        // pretty-print
        xmlKeepBlanksDefault(0);

        // The following two lines should set indenting to 4 spaces but it does
        // not work. 2 spaces are default.
        //
        // xmlThrDefIndentTreeOutput(1);
        // xmlThrDefTreeIndentString("");
        auto doc =
            xmlParseMemory(prj_stream.str().c_str(), prj_stream.str().size());
        auto prj_out =
            std::filesystem::path(prj_file).parent_path() / "_project.prj";
        xmlSaveFormatFileEnc(prj_out.c_str(), doc, "utf-8", 1);
        xmlFreeDoc(doc);
    }
    xmlCleanupParser();
}
}  // namespace BaseLib
