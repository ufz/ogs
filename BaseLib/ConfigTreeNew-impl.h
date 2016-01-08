/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConfigTreeNew.h"

#include <logog/include/logog.hpp>

namespace BaseLib
{

//! Wraps a pair of iterators for use as a range in range-based for-loops.
template<typename Iterator>
class Range
{
public:
    explicit Range(Iterator begin, Iterator end)
        : _begin(begin), _end(end)
    {}

    Iterator begin() const { return _begin; }
    Iterator end()   const { return _end; }
private:
    Iterator _begin;
    Iterator _end;
};

template<typename T>
T
ConfigTreeNew::
getConfParam(std::string const& param) const
{
    auto p = getConfParamOptional<T>(param);
    if (p) return *p;

    error("Key <" + param + "> has not been found");
    return T();
}

template<typename T>
T
ConfigTreeNew::
getConfParam(std::string const& param, T const& default_value) const
{
    auto p = getConfParamOptional<T>(param);
    if (p) return *p;
    return default_value;
}

template<typename T>
boost::optional<T>
ConfigTreeNew::
getConfParamOptional(std::string const& param) const
{
    checkUnique(param);
    auto p = _tree->get_child_optional(param);

    bool peek_only = p == boost::none;
    markVisited<T>(param, peek_only);

    if (p) {
        auto v = p->get_value_optional<T>();
        if (v) {
            return v;
        } else {
            error("Value for key <" + param + "> `" + shortString(p->data())
                  + "' not convertible to the desired type.");
        }
    }

    return boost::none;
}

template<typename T>
Range<ConfigTreeNew::ValueIterator<T> >
ConfigTreeNew::
getConfParamList(std::string const& param) const
{
    checkUnique(param);
    markVisited<T>(param, true);

    auto p = _tree->equal_range(param);
    return Range<ValueIterator<T> >(
                ValueIterator<T>(p.first,  param, *this),
                ValueIterator<T>(p.second, param, *this));
}

template<typename T>
T
ConfigTreeNew::
peekConfParam(std::string const& param) const
{
    checkKeyname(param);

    auto p =_tree->get_child_optional(param);

    if (!p) {
        error("Key <" + param + "> has not been found");
    } else {
        try {
            return p->get_value<T>();
        } catch (boost::property_tree::ptree_bad_data) {
            error("Value for key <" + param + "> `" + shortString(p->data())
                  + "' not convertible to the desired type.");
        }
    }

    return T();
}

template<typename T>
void
ConfigTreeNew::
checkConfParam(std::string const& param, T const& value) const
{
    if (getConfParam<T>(param) != value) {
        error("The value of key <" + param + "> is not the expected one.");
    }
}

template<typename Ch>
void
ConfigTreeNew::
checkConfParam(std::string const& param, Ch const* value) const
{
    if (getConfParam<std::string>(param) != value) {
        error("The value of key <" + param + "> is not the expected one.");
    }
}


template<typename T>
ConfigTreeNew::CountType&
ConfigTreeNew::
markVisited(std::string const& key, bool peek_only) const
{
    auto const type = std::type_index(typeid(T));

    auto p = _visited_params.emplace(key, CountType{peek_only ? 0 : 1, type});

    if (!p.second) { // no insertion happened
        auto& v = p.first->second;
        if (v.type == type) {
            if (!peek_only) ++v.count;
        } else {
            error("There already was an attempt to obtain key <" + key
                  + "> with type \"" + v.type.name() + "\" (now: \"" + type.name() + "\").");
        }
    }

    return p.first->second;
}

}
