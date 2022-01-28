/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ConfigTree.h"

namespace BaseLib
{
/*! Create a ConfigTree from an XML file.
 *
 * \param filepath     see ConfigTreeTopLevel::ConfigTreeTopLevel()
 * \param be_ruthless  see ConfigTreeTopLevel::ConfigTreeTopLevel()
 * \param toplevel_tag name of the outermost tag in the XML file. The returned
 *                     ConfigTree is rooted one level below that tag.
 * \param patch_files  optional vector of strings with patch file paths.
 *
 * The parameter \c toplevel_tag is provided for compatibility with our existing
 * configuration files whose toplevel tags are written in camel case, which
 * conflicts with the naming rules of ConfigTree. Via that parameter the naming
 * rules do not apply to the toplevel tag.
 *
 * Unfortunately the XML parser shipped with <tt>boost::property_tree</tt> does
 * not fully support the XML standard. Additionally there might be encoding
 * issues. From their docs:
 *
 * > Please note that RapidXML does not understand the encoding specification.
 * If you pass it a > character buffer, it assumes the data is already correctly
 * encoded; if you pass it a filename, > it will read the file using the
 * character conversion of the locale you give it (or the global > locale if you
 * give it none). This means that, in order to parse a UTF-8-encoded XML file
 * into > a wptree, you have to supply an alternate locale, either directly or
 * by replacing the global one.
 *
 * \see http://www.boost.org/doc/libs/1_60_0/doc/html/property_tree/parsers.html
 */
ConfigTree makeConfigTree(std::string const& filepath,
                          bool const be_ruthless,
                          std::string const& toplevel_tag,
                          std::stringstream& prj_stream);

ConfigTree makeConfigTreeFromFile(const std::string& filepath,
                                  const bool be_ruthless,
                                  const std::string& toplevel_tag);

}  // namespace BaseLib
