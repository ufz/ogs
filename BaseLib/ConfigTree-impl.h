// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// IWYU pragma: private, include "ConfigTree.h"

#include <exception>
#include <sstream>
#include <utility>
#include <vector>

#include "ConfigTree.h"
#include "DemangleTypeInfo.h"
#include "StringTools.h"

namespace BaseLib
{
//! Wraps a pair of iterators for use as a range in range-based for-loops.
template <typename Iterator>
class Range
{
public:
    explicit Range(Iterator begin, Iterator end)
        : begin_(std::move(begin)), end_(std::move(end))
    {
    }

    Iterator begin() const { return begin_; }
    Iterator end() const { return end_; }
    std::size_t size() const { return std::distance(begin_, end_); }
    bool empty() const { return size() == 0; }

private:
    Iterator begin_;
    Iterator end_;
};

/*! Read-only view of a child produced by getAllChildren().
 *
 * A thin move-only wrapper around the child's own ConfigTree subtree.
 * getAllChildren() marks every child tag as consumed in the parent, so the
 * parent does not additionally warn about the tags themselves; instead each
 * child polices itself.
 *
 * A child counts as read once it is dereferenced via operator->() /
 * operator*(). A child that is never dereferenced warns on destruction as
 * unread whatever its content - this catches a whole child that was silently
 * ignored, e.g. one whose tag name was mistyped. Once dereferenced,
 * the child behaves exactly like a ConfigTree obtained from getConfigSubtree():
 * it owns its subtree and runs the normal checkAndInvalidate() on destruction,
 * so any of its sub-children / attributes / immediate data that were not read
 * produce the usual unread-element warnings.
 *
 * Like a getConfigSubtree() result, a Child keeps the underlying tree alive
 * through its shared top-level pointer, so it may safely outlive the ConfigTree
 * it came from.
 */
class ConfigTree::Child
{
public:
    Child(Child&& other) noexcept
        : ct_(std::move(other.ct_)), accessed_(other.accessed_)
    {
    }
    Child(Child const&) = delete;

    ConfigTree const* operator->() const
    {
        accessed_ = true;
        return &ct_;
    }
    ConfigTree const& operator*() const
    {
        accessed_ = true;
        return ct_;
    }

    /*! A child returned by getAllChildren() but never dereferenced counts as
     * unread and warns on destruction, whatever its content.
     *
     * This is what makes getAllChildren() catch a whole child that was silently
     * ignored - e.g. one whose tag name was mistyped - which the
     * subtree's own content check cannot: a content-less child has nothing to
     * report. A dereferenced child instead defers to that content check (its
     * unread sub-children / attributes / immediate data warn as usual), so the
     * warning here is suppressed and not duplicated.
     */
    ~Child()
    {
        // ct_.tree_ == nullptr means this Child was moved from (the moved-to
        // Child owns the subtree now); std::uncaught_exceptions() > 0 mirrors
        // ConfigTree's own destructor, which stays silent during unwinding.
        if (!accessed_ && ct_.tree_ != nullptr &&
            std::uncaught_exceptions() == 0)
        {
            ct_.warning(
                "This child, returned by getAllChildren(), has not been "
                "read.");
            ct_.tree_ = nullptr;  // skip the subtree's own redundant check
        }
    }

private:
    friend class ConfigTree;
    explicit Child(ConfigTree&& ct) : ct_(std::move(ct)) {}

    ConfigTree ct_;
    //! Whether operator->() / operator*() was ever called. Mutable because both
    //! are const accessors; a Child is logically const while being read.
    mutable bool accessed_ = false;
};

template <typename T>
T ConfigTree::getConfigParameter(std::string const& param) const
{
    if (auto p = getConfigParameterOptional<T>(param))
    {
        return *p;
    }

    error("Key <" + param + "> has not been found");
}

template <typename T>
T ConfigTree::getConfigParameter(std::string const& param,
                                 T const& default_value) const
{
    if (auto p = getConfigParameterOptional<T>(param))
    {
        return *p;
    }

    return default_value;
}

template <typename T>
std::optional<T> ConfigTree::getConfigParameterOptional(
    std::string const& param) const
{
    checkUnique(param);

    return getConfigParameterOptionalImpl(param, static_cast<T*>(nullptr));
}

template <typename T>
std::optional<T> ConfigTree::getConfigParameterOptionalImpl(
    std::string const& param, T* /*unused*/) const
{
    if (auto p = getConfigSubtreeOptional(param))
    {
        return p->getValue<T>();
    }

    return std::nullopt;
}

template <typename T>
std::optional<std::vector<T>> ConfigTree::getConfigParameterOptionalImpl(
    std::string const& param, std::vector<T>* /*unused*/) const
{
    if (auto p = getConfigSubtreeOptional(param))
    {
        std::string const raw = p->getValue<std::string>();

        std::size_t bad_idx = 0;
        if (auto parsed = BaseLib::tryParseVector<T>(raw, &bad_idx))
        {
            return parsed;  // OK
        }

        // preserve original error wording & token index
        error("Value for key <" + param + "> `" + shortString(raw) +
              "' not convertible to a vector of the desired type."
              " Could not convert token no. " +
              std::to_string(bad_idx) + ".");
        return std::nullopt;  // not reached
    }

    return std::nullopt;
}

template <typename T>
Range<ConfigTree::ValueIterator<T>> ConfigTree::getConfigParameterList(
    std::string const& param) const
{
    checkUnique(param);
    markVisited<T>(param, Attr::TAG, true);

    auto p = tree_->equal_range(param);
    return Range<ValueIterator<T>>(ValueIterator<T>(p.first, param, *this),
                                   ValueIterator<T>(p.second, param, *this));
}

template <typename T>
T ConfigTree::peekConfigParameter(std::string const& param) const
{
    checkKeyname(param);

    if (auto p = tree_->get_child_optional(param))
    {
        try
        {
            return p->get_value<T>();
        }
        catch (boost::property_tree::ptree_bad_data const&)
        {
            error("Value for key <" + param + "> `" + shortString(p->data()) +
                  "' not convertible to the desired type.");
        }
    }
    else
    {
        error("Key <" + param + "> has not been found");
    }
}

template <typename T>
T ConfigTree::getValue() const
{
    if (have_read_data_)
    {
        error("The data of this subtree has already been read.");
    }

    have_read_data_ = true;

    if (auto v = tree_->get_value_optional<T>())
    {
        return *v;
    }
    error("Value `" + shortString(tree_->data()) +
          "' is not convertible to the desired type.");
}

template <typename T>
T ConfigTree::getConfigAttribute(std::string const& attr) const
{
    if (auto a = getConfigAttributeOptional<T>(attr))
    {
        return *a;
    }

    error("Did not find XML attribute with name '" + attr + "'.");
}

template <typename T>
T ConfigTree::getConfigAttribute(std::string const& attr,
                                 T const& default_value) const
{
    if (auto a = getConfigAttributeOptional<T>(attr))
    {
        return *a;
    }

    return default_value;
}

template <typename T>
std::optional<T> ConfigTree::getConfigAttributeOptional(
    std::string const& attr) const
{
    checkUniqueAttr(attr);
    auto& ct = markVisited<T>(attr, Attr::ATTR, true);

    if (auto attrs = tree_->get_child_optional("<xmlattr>"))
    {
        if (auto a = attrs->get_child_optional(attr))
        {
            ++ct.count;  // count only if attribute has been found
            if (auto v = a->get_value_optional<T>())
            {
                return std::make_optional(*v);
            }
            error("Value for XML attribute '" + attr + "' `" +
                  shortString(a->data()) +
                  "' not convertible to the desired type.");
        }
    }

    return std::nullopt;
}

template <typename T>
ConfigTree::CountType& ConfigTree::markVisited(std::string const& key,
                                               Attr const is_attr,
                                               bool const peek_only) const
{
    auto const type = std::type_index(typeid(T));

    auto p = visited_params_.emplace(std::make_pair(is_attr, key),
                                     CountType{peek_only ? 0 : 1, type});

    if (!p.second)
    {  // no insertion happened
        auto& v = p.first->second;
        if (v.type == type)
        {
            if (!peek_only)
            {
                ++v.count;
            }
        }
        else
        {
            error("There already was an attempt to obtain key <" + key +
                  "> with type '" + BaseLib::demangle(v.type.name()) +
                  "' (now: '" + type.name() + "').");
        }
    }

    return p.first->second;
}

}  // namespace BaseLib
