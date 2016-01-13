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
	: _tree(other._tree)
	, _path(other._path)
	, _visited_params(std::move(other._visited_params))
	, _onerror(other._onerror)
	, _onwarning(other._onwarning)
{
	other._tree = nullptr;
}

ConfigTreeNew::~ConfigTreeNew()
{
	if (!_tree) return;

	for (auto const& p : *_tree)
	{
		markVisitedDecrement(p.first);
	}

	for (auto const& p : _visited_params)
	{
		if (p.second.count > 0) {
			warning("Key <" + p.first + "> has been read " + std::to_string(p.second.count)
					+ " time(s) more than it was present in the configuration tree.");
		} else if (p.second.count < 0) {
			warning("Key <" + p.first + "> has been read " + std::to_string(-p.second.count)
					+ " time(s) less than it was present in the configuration tree.");
		}
	}
}

ConfigTreeNew
ConfigTreeNew::
getConfSubtree(std::string const& root)
{
    if (auto t = getConfSubtreeOptional(root)) {
        return std::move(*t);
    } else {
        error("Key <" + root + "> has not been found.");
        return ConfigTreeNew(PTree(), *this, ""); // TODO that will crash
    }
}

boost::optional<ConfigTreeNew>
ConfigTreeNew::
getConfSubtreeOptional(std::string const& root)
{
    checkUnique(root);
    auto subtree = _tree->get_child_optional(root);

    if (subtree) {
        markVisited(root);
        return boost::optional<ConfigTreeNew>(std::move(
                ConfigTreeNew(*subtree, *this, root)));
    } else {
        markVisited(root, true);
        return boost::optional<ConfigTreeNew>();
    }
}

Range<ConfigTreeNew::SubtreeIterator>
ConfigTreeNew::
getConfSubtreeList(std::string const& root)
{
    checkUnique(root);
    markVisited(root, true);

    auto p = _tree->equal_range(root);

    return Range<SubtreeIterator>(
                SubtreeIterator(p.first,  root, *this),
                SubtreeIterator(p.second, root, *this));
}

void ConfigTreeNew::ignoreConfParam(const std::string &param)
{
    checkUnique(param);
    // if not found, peek only
    bool peek_only = _tree->find(param) == _tree->not_found();
    markVisited(param, peek_only);
}

void ConfigTreeNew::ignoreConfParamAll(const std::string &param)
{
    checkUnique(param);
    auto& ct = markVisited(param, true);

    auto p = _tree->equal_range(param);
    for (auto it = p.first; it != p.second; ++it) {
        ++ct.count;
    }
}


void ConfigTreeNew::error(const std::string& message)
{
	_onerror(_path, message);
}

void ConfigTreeNew::warning(const std::string& message)
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


void ConfigTreeNew::checkKeyname(std::string const& key)
{
	if (key.empty()) {
		error("Search for empty key.");
	} else if (key_chars_start.find(key.front()) == std::string::npos) {
		error("Key <" + key + "> starts with an illegal character.");
	} else if (key.find_first_not_of(key_chars, 1) != std::string::npos) {
		error("Key <" + key + "> contains illegal characters.");
	}
}

std::string ConfigTreeNew::
joinPaths( const std::string &p1, const std::string &p2)
{
	if (p2.empty()) {
		error("Second path to be joined is empty.");
	}

	if (p1.empty()) return p2;

	return p1 + pathseparator + p2;
}

void ConfigTreeNew::checkUnique(const std::string &key)
{
	checkKeyname(key);

	if (_visited_params.find(key) != _visited_params.end()) {
		error("Key <" + key + "> has already been processed.");
	}
}

ConfigTreeNew::CountType&
ConfigTreeNew::
markVisited(std::string const& key, bool peek_only)
{
    return markVisited<ConfigTreeNew>(key, peek_only);
}

void
ConfigTreeNew::
markVisitedDecrement(std::string const& key)
{
    auto const type = std::type_index(typeid(nullptr));

    auto p = _visited_params.emplace(key, CountType{-1, type});

    if (!p.second) { // no insertion happened
        auto& v = p.first->second;
        --v.count;
    }
}

}
