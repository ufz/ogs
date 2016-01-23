/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConfigTreeNew.h"

namespace BaseLib
{

const char ConfigTreeNew::pathseparator = '/';
const std::string ConfigTreeNew::key_chars_start = "abcdefghijklmnopqrstuvwxyz";
const std::string ConfigTreeNew::key_chars = key_chars_start + "_0123456789";

ConfigTreeNew::
ConfigTreeNew(PTree const& tree,
              Callback const& error_cb,
              Callback const& warning_cb)
    : _tree(&tree), _onerror(error_cb), _onwarning(warning_cb)
{
    if (!_onerror) {
        ERR("ConfigTree: No valid error handler provided.");
        std::abort();
    }
    if (!_onwarning) {
        ERR("ConfigTree: No valid warning handler provided.");
        std::abort();
    }
}

ConfigTreeNew::
ConfigTreeNew(PTree const& tree, ConfigTreeNew const& parent,
              std::string const& root)
    : _tree(&tree), _path(joinPaths(parent._path, root)),
      _onerror(parent._onerror), _onwarning(parent._onwarning)
{
    checkKeyname(root);
}

ConfigTreeNew::
ConfigTreeNew(ConfigTreeNew && other)
    : _tree                    (other._tree)
    , _path          (std::move(other._path))
    , _visited_params(std::move(other._visited_params))
    , _have_read_data          (other._have_read_data)
    , _onerror       (std::move(other._onerror))
    , _onwarning     (std::move(other._onwarning))
{
    other._tree = nullptr;
}

ConfigTreeNew::~ConfigTreeNew()
{
    checkAndInvalidate();
}

ConfigTreeNew&
ConfigTreeNew::
operator=(ConfigTreeNew&& other)
{
    checkAndInvalidate();

    _tree           = other._tree;
    other._tree     = nullptr;
    _path           = std::move(other._path);
    _visited_params = std::move(other._visited_params);
    _have_read_data = other._have_read_data;
    _onerror        = std::move(other._onerror);
    _onwarning      = std::move(other._onwarning);

    return *this;
}

ConfigTreeNew
ConfigTreeNew::
getConfParam(std::string const& root) const
{
    auto ct = getConfSubtree(root);
    if (hasChildren(ct))
        error("Requested parameter <" + root + "> actually is a subtree.");
    return ct;
}

boost::optional<ConfigTreeNew>
ConfigTreeNew::
getConfParamOptional(std::string const& root) const
{
    auto ct = getConfSubtreeOptional(root);
    if (ct && hasChildren(*ct))
        error("Requested parameter <" + root + "> actually is a subtree.");
    return ct;
}

ConfigTreeNew
ConfigTreeNew::
getConfSubtree(std::string const& root) const
{
    if (auto t = getConfSubtreeOptional(root)) {
        return std::move(*t);
    } else {
        error("Key <" + root + "> has not been found.");
    }
}

boost::optional<ConfigTreeNew>
ConfigTreeNew::
getConfSubtreeOptional(std::string const& root) const
{
    checkUnique(root);

    if (auto subtree = _tree->get_child_optional(root)) {
        markVisited(root, false);
        return ConfigTreeNew(*subtree, *this, root);
    } else {
        markVisited(root, true);
        return boost::none;
    }
}

Range<ConfigTreeNew::SubtreeIterator>
ConfigTreeNew::
getConfSubtreeList(std::string const& root) const
{
    checkUnique(root);
    markVisited(root, true);

    auto p = _tree->equal_range(root);

    return Range<SubtreeIterator>(
                SubtreeIterator(p.first,  root, *this),
                SubtreeIterator(p.second, root, *this));
}

void ConfigTreeNew::ignoreConfParam(const std::string &param) const
{
    checkUnique(param);
    // if not found, peek only
    bool peek_only = _tree->find(param) == _tree->not_found();
    markVisited(param, peek_only);
}

void ConfigTreeNew::ignoreConfParamAll(const std::string &param) const
{
    checkUnique(param);
    auto& ct = markVisited(param, true);

    auto p = _tree->equal_range(param);
    for (auto it = p.first; it != p.second; ++it) {
        ++ct.count;
    }
}


void ConfigTreeNew::error(const std::string& message) const
{
    _onerror(_path, message);
    std::abort();
}

void ConfigTreeNew::warning(const std::string& message) const
{
    _onwarning(_path, message);
}


void ConfigTreeNew::onerror(const std::string& path, const std::string& message)
{
    ERR("ConfigTree: At path <%s>: %s", path.c_str(), message.c_str());
    std::abort();
}

void ConfigTreeNew::onwarning(const std::string& path, const std::string& message)
{
    WARN("ConfigTree: At path <%s>: %s", path.c_str(), message.c_str());
}

std::string ConfigTreeNew::shortString(const std::string &s)
{
    const std::size_t maxlen = 100;

    if (s.size() < maxlen) return s;

    return s.substr(0, maxlen-3) + "...";
}


void ConfigTreeNew::checkKeyname(std::string const& key) const
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

std::string ConfigTreeNew::
joinPaths( const std::string &p1, const std::string &p2) const
{
    if (p2.empty()) {
        error("Second path to be joined is empty.");
    }

    if (p1.empty()) return p2;

    return p1 + pathseparator + p2;
}

void ConfigTreeNew::checkUnique(const std::string &key) const
{
    checkKeyname(key);

    if (_visited_params.find({false, key}) != _visited_params.end()) {
        error("Key <" + key + "> has already been processed.");
    }
}

void ConfigTreeNew::checkUniqueAttr(const std::string &attr) const
{
    checkKeyname(attr);

    if (_visited_params.find({true, attr}) != _visited_params.end()) {
        error("Attribute \"" + attr + "\" has already been processed.");
    }
}

ConfigTreeNew::CountType&
ConfigTreeNew::
markVisited(std::string const& key, bool const peek_only) const
{
    return markVisited<ConfigTreeNew>(key, false, peek_only);
}

void
ConfigTreeNew::
markVisitedDecrement(bool const is_attr, std::string const& key) const
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
ConfigTreeNew::hasChildren(ConfigTreeNew const& ct) const
{
    auto const& tree = *ct._tree;
    if (tree.begin() == tree.end())
        return false; // no children
    if (tree.front().first == "<xmlattr>"
        && (++tree.begin()) == tree.end())
        return false; // only attributes

    return true;
}

void
ConfigTreeNew::checkAndInvalidate()
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
            markVisitedDecrement(false, p.first);
    }

    // iterate over attributes
    if (auto attrs = _tree->get_child_optional("<xmlattr>")) {
        for (auto const& p : *attrs) {
            markVisitedDecrement(true, p.first);
        }
    }

    for (auto const& p : _visited_params)
    {
        auto const& tag   = p.first.second;
        auto const& count = p.second.count;

        if (p.first.first) { // tag
            if (count > 0) {
                warning("XML attribute \"" + tag + "\" has been read " + std::to_string(count)
                        + " time(s) more than it was present in the configuration tree.");
            } else if (count < 0) {
                warning("XML attribute \"" + tag + "\" has been read " + std::to_string(-count)
                        + " time(s) less than it was present in the configuration tree.");
            }
        } else { // XML attribute
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


void checkAndInvalidate(ConfigTreeNew &conf)
{
    conf.checkAndInvalidate();
}

void checkAndInvalidate(ConfigTreeNew* const conf)
{
    if (conf) conf->checkAndInvalidate();
}

void checkAndInvalidate(std::unique_ptr<ConfigTreeNew> const& conf)
{
    if (conf) conf->checkAndInvalidate();
}

}
