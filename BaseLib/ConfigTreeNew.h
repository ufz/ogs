/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ConfigTree.h"

#include <typeindex>
#include <map>

#include <functional>

namespace BaseLib
{

template<typename Iterator> class Range;

/*!
 * Wrapper around a Boost Property Tree with some basic error reporting features.
 *
 * Features. This class:
 *  * makes sure that every configuration setting in a Property Tree is read
 *    exactly once. If some settings is not read (e.g. due to a typo), a warning message
 *    is generated. The message contains a hint where it occured.
 *  * enforces a naming scheme of settings: letters a-z, numbers 0-9, underscore
 *  * provides some functionality to read lists of values using range-based for loops.
 *  * has rather long method names that are easily greppable from the source code. So a list
 *    of supported configuration options can be easily obtained from the source code.
 *
 * The purpose of this class is to reduce or completely avoid the amount of error-handling
 * code in routines that take configuration parameters.
 *
 * Most methods of this class check that they have not been called before for the same
 * \c ConfigTree and the same parameter. This behaviour helps to enforce that every parameter
 * is read exactly once during parsing of the configuration settings.
 *
 * The most notable restriction of this class when compared to plain tree traversal is, that
 * one must know all the XML tags (i.e. configuration parameters) at compile time. It is not
 * possible to read from this class, which configuration parameters are present in the tree.
 * This restriction, however, is intended, because it provides the possibility to get all
 * existing configuration parameters from the source code.
 */
class ConfigTreeNew final
{
public:
    /*!
     * A wrapper around a Boost Iterator for iterating over ranges of subtrees.
     *
     * The methods of this class tell the associated (parent) \c ConfigTree object when
     * a setting has been parsed.
     */
    class SubtreeIterator
            : public std::iterator<std::input_iterator_tag, ConfigTreeNew>
    {
    public:
        using Iterator = boost::property_tree::ptree::const_assoc_iterator;

        explicit SubtreeIterator(Iterator it, std::string const& root,
                                 ConfigTreeNew& parent)
            : _it(it), _root(root), _parent(parent)
        {}

        SubtreeIterator& operator++() {
            ++_it;
            _has_incremented = true;
            return *this;
        }

        ConfigTreeNew operator*() {
            // if this iterator has been incremented since the last dereference,
            // tell the _parent instance that a subtree now has been parsed.
            if (_has_incremented) {
                _has_incremented = false;
                _parent.markVisited(_root);
            }
            return ConfigTreeNew(_it->second, _parent, _root);
        }

        bool operator==(SubtreeIterator const& other) const {
            return _it == other._it;
        }

        bool operator!=(SubtreeIterator const& other) const {
            return _it != other._it;
        }

    private:
        bool _has_incremented = true;
        Iterator _it;
        std::string const _root;
        ConfigTreeNew& _parent;
    };


    /*!
     * A wrapper around a Boost Iterator for iterating over ranges of values.
     *
     * The methods of this class tell the associated (parent) \c ConfigTree object when
     * a setting has been parsed.
     */
    template<typename ValueType>
    class ValueIterator
            : public std::iterator<std::input_iterator_tag, ValueType>
    {
    public:
        using Iterator = boost::property_tree::ptree::const_assoc_iterator;

        explicit ValueIterator(Iterator it, std::string const& root,
                               ConfigTreeNew& parent)
            : _it(it), _root(root), _parent(parent)
        {}

        ValueIterator<ValueType>& operator++() {
            ++_it;
            _has_incremented = true;
            return *this;
        }

        ValueType operator*() {
            // if this iterator has been incremented since the last dereference,
            // tell the _parent instance that a setting now has been parsed.
            if (_has_incremented) {
                _has_incremented = false;
                _parent.markVisited<ValueType>(_root);
            }

            if (_it->second.begin() != _it->second.end()) {
                _parent.error("Configuration at key " + _root + " has subitems.");
                return ValueType();
            }

            auto v = _it->second.get_value_optional<ValueType>();

            if (v) return *v;

            // TODO: change error method
            _parent.error("Could not get value out of key " + _root + ".");
            return ValueType();
        }

        bool operator==(ValueIterator<ValueType> const& other) const {
            return _it == other._it;
        }

        bool operator!=(ValueIterator<ValueType> const& other) const {
            return _it != other._it;
        }

    private:
        bool _has_incremented = true;
        Iterator _it;
        std::string const _root;
        ConfigTreeNew& _parent;
    };

    using PTree = boost::property_tree::ptree;

    //! Type of the function objects used as callbacks.
    //! The first argument denotes the path in the tree at which an event (warning/error)
    //! occured, the second argument is the associated message
    using Callback = std::function<void(const std::string& path,
                                        const std::string& message)>;

    /*! Creates a new instance wrapping the given Boost Property Tree.
     *
     * \param tree the Boost Property Tree to be wrapped
     * \param error_cb callback function to be called on error.
     * \param warning_cb callback function to be called on warning.
     *
     * The callback functions must be valid callable functions, i.e. not nullptr's.
     * They are configurable in order to make unit tests of this class easier.
     * They should not be provided in production code!
     *
     * If a custom error callback is provided, this function should break out of
     * the normal execution order, e.g., by throwing or by calling std::abort(),
     * because otherwise this class will effectively treat errors as no-errors.
     */
    explicit ConfigTreeNew(PTree const& tree,
                           Callback const& error_cb = onerror,
                           Callback const& warning_cb = onwarning);

    //! copying is not compatible with the semantics of this class
    ConfigTreeNew(ConfigTreeNew const&) = delete;

    //! After being moved from, \c other is in an undefined state and must not be
    //! used anymore!
    ConfigTreeNew(ConfigTreeNew && other);

    ConfigTreeNew() = delete;

    void operator=(ConfigTreeNew const&) = delete;
    void operator=(ConfigTreeNew &&) = delete;

    /*! Get parameter \c param of type \c T from the configuration tree.
     *
     * \return the value looked for or a default constructed value \c T() in the case of
     *         an error. For the behaviour in case of an error, see also the documentation
     *         of the method error().
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    template<typename T> T
    getConfParam(std::string const& param);

    /*! Get parameter \c param of type \c T from the configuration tree or the \c default_value.
     *
     * This method has a similar behaviour as getConfParam(std::string const&) except in case
     * of errors the \c default_value is returned.
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    template<typename T> T
    getConfParam(std::string const& param, T const& default_value);

    /*! Get parameter \c param of type \c T from the configuration tree if present
     *
     * This method has a similar behaviour as getConfParam(std::string const&) except
     * no errors are raised. Rather it can be told from the return value if the
     * parameter could be read.
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    template<typename T> boost::optional<T>
    getConfParamOptional(std::string const& param);

    /*! Returns all parameters with the name \c param from the current level of the tree.
     *
     * The return value is suitable to be used with range-base for-loops.
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    template<typename T> Range<ValueIterator<T> >
    getConfParamList(std::string const& param);

    /*! Peek at a parameter \c param of type \c T from the configuration tree.
     *
     * This method is an exception to the single-read rule. It is meant to be used to
     * tell from a ConfigTree instance where to pass that instance on for further processing.
     *
     * Return value and error behaviour are the same as for getConfParam(std::string const&).
     */
    template<typename T> T
    peekConfParam(std::string const& param);

    /*! Assert that \c param has the given \c value.
     *
     * Convenience method combining getConfParam(std::string const&) with a check.
     */
    template<typename T> void
    checkConfParam(std::string const& param, T const& value);

    //! Make checkConfParam() work for string literals.
    template<typename Ch> void
    checkConfParam(std::string const& param, Ch const* value);

    /*! Get the subtree rooted at \c root
     *
     * If \c root is not found error() is called.
     *
     * \pre \c root must not have been read before from this ConfigTree.
     */
    ConfigTreeNew
    getConfSubtree(std::string const& root);

    /*! Get the subtree rooted at \c root if present
     *
     * \pre \c root must not have been read before from this ConfigTree.
     */
    boost::optional<ConfigTreeNew>
    getConfSubtreeOptional(std::string const& root);

    /*! Get all subtrees that have a root \c root from the current level of the tree.
     *
     * The return value is suitable to be used with range-base for-loops.
     *
     * \pre \c root must not have been read before from this ConfigTree.
     */
    Range<SubtreeIterator>
    getConfSubtreeList(std::string const& root);

    /*! Tell this instance to ignore parameter \c param.
     *
     * This method is used to avoid warning messages.
     *
     * \pre \c root must not have been read before from this ConfigTree.
     */
    void ignoreConfParam(std::string const& param);

    /*! Tell this instance to ignore all parameters \c param on the current level of the tree.
     *
     * This method is used to avoid warning messages.
     *
     * \pre \c root must not have been read before from this ConfigTree.
     */
    void ignoreConfParamAll(std::string const& param);

    //! The destructor performs the check if all nodes at the current level of the tree
    //! have been read.
    ~ConfigTreeNew();

private:
    struct CountType
    {
        int count;
        std::type_index type;
    };

    //! Used for wrapping a subtree
    explicit ConfigTreeNew(PTree const& tree, ConfigTreeNew const& parent, std::string const& root);

    //! Called if an error occurs. Will call the error callback.
    //! This method only acts as a helper method.
    void error(std::string const& message);

    //! Called for printing warning messages. Will call the warning callback.
    //! This method only acts as a helper method.
    void warning(std::string const& message);

    //! Checks if \c key complies with the rules [a-z0-9_].
    void checkKeyname(std::string const& key);

    //! Used to generate the path of a subtree.
    std::string joinPaths(std::string const& p1, std::string const& p2);

    //! Asserts that the \c key has not been read yet.
    void checkUnique(std::string const& key);

    /*! Keeps track of the key \c key and its value type \c T.
     *
     * This method asserts that a key is read always with the same type.
     *
     * \c param peek_only if true, do not change the read-count of the given key.
     */
    template<typename T>
    CountType& markVisited(std::string const& key, bool peek_only = false);

    /*! Keeps track of the key \c key and its value type ConfigTree.
     *
     * This method asserts that a key is read always with the same type.
     *
     * \c param peek_only if true, do not change the read-count of the given key.
     */
    CountType& markVisited(std::string const& key, bool peek_only = false);

    //! Used in the destructor to compute the difference between number of reads of a parameter
    //! and the number of times it exists in the ConfigTree
    void markVisitedDecrement(std::string const& key);

    //! Default error callback function
    //! Will print an error message and call std::abort()
    static void onerror(std::string const& path, std::string const& message);

    //! Default warning callback function
    //! Will print a warning message
    static void onwarning(std::string const& path, std::string const& message);

    //! returns a short string at suitable for error/warning messages
    static std::string shortString(std::string const& s);

    //! The wrapped tree.
    boost::property_tree::ptree const* _tree;

    //! A path printed in error/warning messages.
    std::string const _path;

    //! A map key -> (count, type) keeping track which parameters have been read how often
    //! and which datatype they have.
    std::map<std::string, CountType> _visited_params;

    const Callback _onerror;
    const Callback _onwarning;

    //! Character separating two path components.
    static const char pathseparator;

    //! Set of allowed characters as the first letter of a key name.
    static const std::string key_chars_start;

    //! Set of allowed characters in a key name.
    static const std::string key_chars;
};

}

#include "ConfigTreeNew-impl.h"

