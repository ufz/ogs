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
    if (auto p = getConfParamOptional<T>(param))
        return *p;

    error("Key <" + param + "> has not been found");
}

template<typename T>
T
ConfigTreeNew::
getConfParam(std::string const& param, T const& default_value) const
{
    if (auto p = getConfParamOptional<T>(param))
        return *p;

    return default_value;
}

template<typename T>
boost::optional<T>
ConfigTreeNew::
getConfParamOptional(std::string const& param) const
{
    checkUnique(param);

    if (auto p = getConfSubtreeOptional(param))
        return p->getValue<T>();

    return boost::none;
}

template<typename T>
Range<ConfigTreeNew::ValueIterator<T> >
ConfigTreeNew::
getConfParamList(std::string const& param) const
{
    checkUnique(param);
    markVisited<T>(param, Attr::TAG, true);

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

    if (auto p =_tree->get_child_optional(param)) {
        try {
            return p->get_value<T>();
        } catch (boost::property_tree::ptree_bad_data) {
            error("Value for key <" + param + "> `" + shortString(p->data())
                  + "' not convertible to the desired type.");
        }
    } else {
        error("Key <" + param + "> has not been found");
    }
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
T
ConfigTreeNew::
getValue() const
{
    if (_have_read_data) {
        error("The data of this subtree has already been read.");
    }

    _have_read_data = true;

    if (auto v = _tree->get_value_optional<T>()) {
        return *v;
    } else {
        error("Value `" + shortString(_tree->data())
              + "' is not convertible to the desired type.");
    }
}

template<typename T>
T
ConfigTreeNew::
getConfAttribute(std::string const& attr) const
{
    if (auto a = getConfAttributeOptional<T>(attr))
        return *a;

    error("Did not find XML attribute with name \"" + attr + "\".");
}

template<typename T>
boost::optional<T>
ConfigTreeNew::
getConfAttributeOptional(std::string const& attr) const
{
    checkUniqueAttr(attr);
    auto& ct = markVisited<T>(attr, Attr::ATTR, true);

    if (auto attrs = _tree->get_child_optional("<xmlattr>")) {
        if (auto a = attrs->get_child_optional(attr)) {
            ++ct.count; // count only if attribute has been found
            if (auto v = a->get_value_optional<T>()) {
                return v;
            } else {
                error("Value for XML attribute \"" + attr + "\" `"
                      + shortString(a->data())
                      + "' not convertible to the desired type.");
            }
        }
    }

    return boost::none;
}

template<typename T>
ConfigTreeNew::CountType&
ConfigTreeNew::
markVisited(std::string const& key, Attr const is_attr,
            bool const peek_only) const
{
    auto const type = std::type_index(typeid(T));

    auto p = _visited_params.emplace(std::make_pair(is_attr, key),
                                     CountType{peek_only ? 0 : 1, type});

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
