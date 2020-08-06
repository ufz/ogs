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
#include <regex>

#include "Error.h"
#include "filesystem.h"
#include "Logging.h"

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

// From https://stackoverflow.com/a/1567703/80480
class line {
    std::string data;
public:
    friend std::istream &operator>>(std::istream &is, line &l) {
        std::getline(is, l.data);
        return is;
    }
    operator std::string() const { return data; }
};

ConfigTreeTopLevel
makeConfigTree(const std::string& filepath, const bool be_ruthless,
               const std::string& toplevel_tag)
{
    ConfigTree::PTree ptree;

    // note: Trimming whitespace and ignoring comments is crucial in order
    //       for our configuration tree implementation to work!
    try
    {
        // Searching for <include filename=".." /> - tags and replacing by the
        // content of the included file.
        std::ifstream file (filepath);
        std::stringstream output;

        std::transform( std::istream_iterator<line>( file ),
            std::istream_iterator<line>(),
            std::ostream_iterator<std::string>( output, "\n" ),
            [](std::string const& _line)
            {
                const std::regex base_regex(".*<include file=\"(.*)\" ?/>.*");
                std::smatch base_match;
                if (std::regex_match(_line, base_match, base_regex))
                {
                    if (base_match.size() == 2)
                    {
                        std::string include_filename = base_match[1].str();
                        std::ifstream include_file (include_filename);
                        std::string include_content(
                            (std::istreambuf_iterator<char>(include_file) ),
                            (std::istreambuf_iterator<char>()) );
                        INFO("Including {:s} into project file.",
                             include_filename);
                        return include_content;
                    }
                }
                return _line;
            }
        );

        read_xml(output, ptree,
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
