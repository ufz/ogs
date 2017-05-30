/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConfigTree.h"

#include <forward_list>
#include <logog/include/logog.hpp>
#include <utility>

#include "Error.h"

// Explicitly instantiate the boost::property_tree::ptree which is a typedef to
// the following basic_ptree.
template class boost::property_tree::basic_ptree<std::string, std::string,
                                                 std::less<std::string>>;

//! Collects swallowed error messages raised by the check during destruction of
//! ConfigTree instances.
static std::forward_list<std::string> configtree_destructor_error_messages;

namespace BaseLib
{

const char ConfigTree::pathseparator = '/';
const std::string ConfigTree::key_chars_start = "abcdefghijklmnopqrstuvwxyz";
const std::string ConfigTree::key_chars = key_chars_start + "_0123456789";

ConfigTree::ConfigTree(PTree const& tree,
                       std::string filename,
                       Callback error_cb,
                       Callback warning_cb)
    : _tree(&tree),
      _filename(std::move(filename)),
      _onerror(std::move(error_cb)),
      _onwarning(std::move(warning_cb))
{
    if (!_onerror) {
        OGS_FATAL("ConfigTree: No valid error handler provided.");
    }
    if (!_onwarning) {
        OGS_FATAL("ConfigTree: No valid warning handler provided.");
    }
}

ConfigTree::
ConfigTree(PTree const& tree, ConfigTree const& parent,
              std::string const& root)
    : _tree(&tree), _path(joinPaths(parent._path, root)),
      _filename(parent._filename),
      _onerror(parent._onerror), _onwarning(parent._onwarning)
{
    checkKeyname(root);
}

ConfigTree::
ConfigTree(ConfigTree && other)
    : _tree                    (other._tree)
    , _path          (std::move(other._path))
    , _filename      (std::move(other._filename))
    , _visited_params(std::move(other._visited_params))
    , _have_read_data          (other._have_read_data)
    , _onerror       (std::move(other._onerror))
    , _onwarning     (std::move(other._onwarning))
{
    other._tree = nullptr;
}

ConfigTree::~ConfigTree()
{
    if (std::uncaught_exception()) {
        /* If the stack unwinds the check below shall be suppressed in order to
         * not accumulate false-positive configuration errors.
         */
        return;
    }

    try {
        checkAndInvalidate();
    } catch (std::exception& e) {
        ERR("%s", e.what());
        configtree_destructor_error_messages.push_front(e.what());
    }
}

ConfigTree&
ConfigTree::
operator=(ConfigTree&& other)
{
    checkAndInvalidate();

    _tree           = other._tree;
    other._tree     = nullptr;
    _path           = std::move(other._path);
    _filename       = std::move(other._filename);
    _visited_params = std::move(other._visited_params);
    _have_read_data = other._have_read_data;
    _onerror        = std::move(other._onerror);
    _onwarning      = std::move(other._onwarning);

    return *this;
}

ConfigTree
ConfigTree::
getConfigParameter(std::string const& root) const
{
    auto ct = getConfigSubtree(root);
    if (ct.hasChildren())
        error("Requested parameter <" + root + "> actually is a subtree.");
    return ct;
}

boost::optional<ConfigTree>
ConfigTree::
getConfigParameterOptional(std::string const& root) const
{
    auto ct = getConfigSubtreeOptional(root);
    if (ct && ct->hasChildren())
        error("Requested parameter <" + root + "> actually is a subtree.");
    return ct;
}

Range<ConfigTree::ParameterIterator>
ConfigTree::
getConfigParameterList(const std::string &param) const
{
    checkUnique(param);
    markVisited(param, Attr::TAG, true);

    auto p = _tree->equal_range(param);

    return Range<ParameterIterator>(
                ParameterIterator(p.first,  param, *this),
                ParameterIterator(p.second, param, *this));
}

ConfigTree
ConfigTree::
getConfigSubtree(std::string const& root) const
{
    if (auto t = getConfigSubtreeOptional(root)) {
        return std::move(*t);
    }
    error("Key <" + root + "> has not been found.");
}

boost::optional<ConfigTree>
ConfigTree::
getConfigSubtreeOptional(std::string const& root) const
{
    checkUnique(root);

    if (auto subtree = _tree->get_child_optional(root)) {
        markVisited(root, Attr::TAG, false);
        return ConfigTree(*subtree, *this, root);
    }
    markVisited(root, Attr::TAG, true);
    return boost::none;
}

Range<ConfigTree::SubtreeIterator>
ConfigTree::
getConfigSubtreeList(std::string const& root) const
{
    checkUnique(root);
    markVisited(root, Attr::TAG, true);

    auto p = _tree->equal_range(root);

    return Range<SubtreeIterator>(
                SubtreeIterator(p.first,  root, *this),
                SubtreeIterator(p.second, root, *this));
}

void ConfigTree::ignoreConfigParameter(const std::string &param) const
{
    checkUnique(param);
    // if not found, peek only
    bool peek_only = _tree->find(param) == _tree->not_found();
    markVisited(param, Attr::TAG, peek_only);
}

void ConfigTree::ignoreConfigAttribute(const std::string &attr) const
{
    checkUniqueAttr(attr);

    // Exercise: Guess what not! (hint: if not found, peek only)
    // Btw. (not a hint) _tree->find() does not seem to work here.
    bool peek_only = !_tree->get_child_optional("<xmlattr>." + attr);

    markVisited(attr, Attr::ATTR, peek_only);
}

void ConfigTree::ignoreConfigParameterAll(const std::string &param) const
{
    checkUnique(param);
    auto& ct = markVisited(param, Attr::TAG, true);

    auto p = _tree->equal_range(param);
    for (auto it = p.first; it != p.second; ++it) {
        ++ct.count;
    }
}


void ConfigTree::error(const std::string& message) const
{
    _onerror(_filename, _path, message);
    OGS_FATAL(
        "ConfigTree: The error handler does not break out of the normal "
        "control flow.");
}

void ConfigTree::warning(const std::string& message) const
{
    _onwarning(_filename, _path, message);
}


void ConfigTree::onerror(const std::string& filename, const std::string& path,
                            const std::string& message)
{
    OGS_FATAL("ConfigTree: In file `%s' at path <%s>: %s",
        filename.c_str(), path.c_str(), message.c_str());
}

void ConfigTree::onwarning(const std::string& filename, const std::string& path,
                              const std::string& message)
{
    WARN("ConfigTree: In file `%s' at path <%s>: %s",
         filename.c_str(), path.c_str(), message.c_str());
}

void ConfigTree::assertNoSwallowedErrors()
{
    if (configtree_destructor_error_messages.empty())
        return;

    ERR("ConfigTree: There have been errors when parsing the configuration "
        "file(s):");

    for (auto const& msg : configtree_destructor_error_messages) {
        ERR("%s", msg.c_str());
    }

    OGS_FATAL("There have been errors when parsing the configuration file(s).");
}

std::string ConfigTree::shortString(const std::string &s)
{
    const std::size_t maxlen = 100;

    if (s.size() < maxlen) return s;

    return s.substr(0, maxlen-3) + "...";
}


void ConfigTree::checkKeyname(std::string const& key) const
{
    if (key.empty()) {
        error("Search for empty key.");
    } else if (key_chars_start.find(key.front()) == std::string::npos) {
        error("Key <" + key + "> starts with an illegal character.");
    } else if (key.find_first_not_of(key_chars, 1) != std::string::npos) {
        error("Key <" + key + "> contains illegal characters.");
    } else if (key.find("__") != std::string::npos) {
        // This is illegal because we use parameter names to generate doxygen
        // page names. Thereby "__" acts as a separator character. Choosing
        // other separators is not possible because of observed limitations
        // for valid doxygen page names.
        error("Key <" + key + "> contains double underscore.");
    }
}

std::string ConfigTree::
joinPaths( const std::string &p1, const std::string &p2) const
{
    if (p2.empty()) {
        error("Second path to be joined is empty.");
    }

    if (p1.empty()) return p2;

    return p1 + pathseparator + p2;
}

void ConfigTree::checkUnique(const std::string &key) const
{
    checkKeyname(key);

    if (_visited_params.find({Attr::TAG, key}) != _visited_params.end()) {
        error("Key <" + key + "> has already been processed.");
    }
}

void ConfigTree::checkUniqueAttr(const std::string &attr) const
{
    // Workaround for handling attributes with xml namespaces and uppercase letters.
    if (attr.find(':') != attr.npos) {
        auto pos = decltype(attr.npos){0};

        // Replace colon and uppercase letters with an allowed character 'a'.
        // That means, attributes containing a colon are also allowed to contain
        // uppercase letters.
        auto attr2 = attr;
        do {
            pos = attr2.find_first_of(":ABCDEFGHIJKLMNOPQRSTUVWXYZ", pos);
            if (pos != attr.npos) attr2[pos] = 'a';
        } while (pos != attr.npos);

        checkKeyname(attr2);
    } else {
        checkKeyname(attr);
    }

    if (_visited_params.find({Attr::ATTR, attr}) != _visited_params.end()) {
        error("Attribute \"" + attr + "\" has already been processed.");
    }
}

ConfigTree::CountType&
ConfigTree::
markVisited(std::string const& key, Attr const is_attr, bool const peek_only) const
{
    return markVisited<ConfigTree>(key, is_attr, peek_only);
}

void
ConfigTree::
markVisitedDecrement(Attr const is_attr, std::string const& key) const
{
    auto const type = std::type_index(typeid(nullptr));

    auto p = _visited_params.emplace(std::make_pair(is_attr, key),
                                     CountType{-1, type});

    if (!p.second) { // no insertion happened
        auto& v = p.first->second;
        --v.count;
    }
}

bool
ConfigTree::hasChildren() const
{
    auto const& tree = *_tree;
    if (tree.begin() == tree.end())
        return false; // no children
    if (tree.front().first == "<xmlattr>"
        && (++tree.begin()) == tree.end())
        return false; // only attributes

    return true;
}

void
ConfigTree::checkAndInvalidate()
{
    if (!_tree) return;

    // Note: due to a limitation in boost::property_tree it is not possible
    // to discriminate between <tag></tag> and <tag/> in the input file.
    // In both cases data() will be empty.
    if ((!_have_read_data) && !_tree->data().empty()) {
        warning("The immediate data `" + shortString(_tree->data())
                +"' of this tag has not been read.");
    }

    // iterate over children
    for (auto const& p : *_tree) {
        if (p.first != "<xmlattr>") // attributes are handled below
            markVisitedDecrement(Attr::TAG, p.first);
    }

    // iterate over attributes
    if (auto attrs = _tree->get_child_optional("<xmlattr>")) {
        for (auto const& p : *attrs) {
            markVisitedDecrement(Attr::ATTR, p.first);
        }
    }

    for (auto const& p : _visited_params)
    {
        auto const& tag   = p.first.second;
        auto const& count = p.second.count;

        switch (p.first.first) {
        case Attr::ATTR:
            if (count > 0) {
                warning("XML attribute \"" + tag + "\" has been read " + std::to_string(count)
                        + " time(s) more than it was present in the configuration tree.");
            } else if (count < 0) {
                warning("XML attribute \"" + tag + "\" has been read " + std::to_string(-count)
                        + " time(s) less than it was present in the configuration tree.");
            }
            break;
        case Attr::TAG:
            if (count > 0) {
                warning("Key <" + tag + "> has been read " + std::to_string(count)
                        + " time(s) more than it was present in the configuration tree.");
            } else if (count < 0) {
                warning("Key <" + tag + "> has been read " + std::to_string(-count)
                        + " time(s) less than it was present in the configuration tree.");
            }
        }
    }

    // The following invalidates this instance, s.t. it can not be read from it anymore,
    // but it also prevents double-checking.
    _tree = nullptr;
}


void checkAndInvalidate(ConfigTree &conf)
{
    conf.checkAndInvalidate();
}

void checkAndInvalidate(ConfigTree* const conf)
{
    if (conf) conf->checkAndInvalidate();
}

void checkAndInvalidate(std::unique_ptr<ConfigTree> const& conf)
{
    if (conf) conf->checkAndInvalidate();
}

}
