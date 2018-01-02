/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ConfigTree.h"

namespace BaseLib
{

/*! Manages a ConfigTree and the <tt>boost::property_tree</tt> it depends on.
 *
 * The whole purpose of this class is making the management of said dependency easy.
 */
class ConfigTreeTopLevel final
{
public:
    /*! Construct a new instance from the given data.
     *
     * \param filepath stored for use in error/warning messages
     * \param be_ruthless if true, then warnings will raise errors, .i.e. lead to program abortion,
     *                    else warnings will only warn
     * \param ptree the underlying ptree of the created ConfigTree
     */
    explicit ConfigTreeTopLevel(std::string const& filepath,
                                bool const be_ruthless,
                                ConfigTree::PTree&& ptree);

    /*! Access the contained ConfigTree.
     *
     * The non-const version of this method has not been implemented in order to prevent invalidating
     * the \c _ctree when it is passed around. In order to check and invalidate \c _ctree use the provided
     * member function.
     */
    ConfigTree const& operator*() const;

    /*! Access the contained ConfigTree.
     *
     * The non-const version of this method has not been implemented in order to prevent invalidating
     * the \c _ctree when it is passed around. In order to check and invalidate \c _ctree use the provided
     * member function.
     */
    ConfigTree const* operator->() const;

    /*! Check if the contained ConfigTree has been processed entirely.
     *
     * This only checks the top level, as usual with ConfigTree instances.
     *
     * \post Afterwards the contained ConfigTree instance must not be used anymore!
     */
    void checkAndInvalidate();

private:
    ConfigTree::PTree const _ptree; //!< <tt>boost::property_tree</tt> that underlies \c _ctree
    ConfigTree              _ctree; //!< ConfigTree depending on \c _ptree
};

/*! Create a ConfigTree from an XML file.
 *
 * \param filepath     see ConfigTreeTopLevel::ConfigTreeTopLevel()
 * \param be_ruthless  see ConfigTreeTopLevel::ConfigTreeTopLevel()
 * \param toplevel_tag name of the outermost tag in the XML file. The returned ConfigTree is rooted
 *                     one level below that tag.
 *
 * The parameter \c toplevel_tag is provided for compatibility with our existing configuration
 * files whose toplevel tags are written in camel case, which conflicts with the naming rules of
 * ConfigTree. Via that parameter the naming rules do not apply to the toplevel tag.
 *
 * Unfortunately the XML parser shipped with <tt>boost::property_tree</tt> does not fully support
 * the XML standard. Additionally there might be encoding issues. From their docs:
 *
 * > Please note that RapidXML does not understand the encoding specification. If you pass it a
 * > character buffer, it assumes the data is already correctly encoded; if you pass it a filename,
 * > it will read the file using the character conversion of the locale you give it (or the global
 * > locale if you give it none). This means that, in order to parse a UTF-8-encoded XML file into
 * > a wptree, you have to supply an alternate locale, either directly or by replacing the global one.
 *
 * \see http://www.boost.org/doc/libs/1_60_0/doc/html/property_tree/parsers.html
 */
ConfigTreeTopLevel
makeConfigTree(std::string const& filepath, bool const be_ruthless,
               std::string const& toplevel_tag);

}
