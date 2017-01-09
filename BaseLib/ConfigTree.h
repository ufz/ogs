/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <typeindex>
#include <map>

#include <functional>
#include <memory>
#include <vector>

#include <boost/property_tree/ptree.hpp>

extern template class boost::property_tree::basic_ptree<
    std::string, std::string, std::less<std::string>>;


namespace BaseLib
{

class ConfigTree;

/*! Check if \c conf has been read entirely and invalidate it.
 *
 * This method can savely be called on \c nullptr's.
 *
 * \see ConfigTree::checkAndInvalidate()
 */
void checkAndInvalidate(ConfigTree* const conf);

//! \overload
void checkAndInvalidate(std::unique_ptr<ConfigTree> const& conf);

//! \overload
void checkAndInvalidate(ConfigTree& conf);

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
 *
 * This class maintains a read counter for each parameter accessed through any of its methods.
 * Read counters are increased with every read (the only exception being the peekConfigParameter() method).
 * The destructor finally decreases the read counter for every tag/attribute it find on the
 * current level of the XML tree. If the increases/decreases don't cancel each other, warning
 * messages are generated. This check can also be enforced before destruction by using the
 * BaseLib::checkAndInvalidate() functions.
 *
 * The design of this class entails some limitations compared to traversing a plain tree,
 * e.g., it is not possible to obtain a list of tags or attributes from the tree,
 * but one has to explicitly query the specific tags/attributes one is interested in.
 * That way it is possible to get all used configuration parameters directly from the source
 * code where this class is used, and to maintain the quality of the configuration parameter
 * documentation.
 *
 * Instances of this class only keep a reference to the underlying <tt>boost::property_tree</tt>.
 * Therefore it is necessary that the underlying property tree stays intact as long as any
 * instance---i.e. the top level ConfigTree and any of its children---reference it.
 * In order to simplify the handling of this dependence, the class ConfigTreeTopLevel can be used.
 *
 * The construction of a ConfigTree from the content of an XML file can be done with the
 * function BaseLib::makeConfigTree(), which performs many error checks. For limitations
 * of the used XML parser, please have a look at that function's documentation.
 */
class ConfigTree final
{
public:
    /*! A wrapper around a Boost Iterator for iterating over ranges of subtrees.
     *
     * The methods of this class tell the associated (parent) \c ConfigTree object when
     * a setting has been parsed.
     */
    class SubtreeIterator
            : public std::iterator<std::input_iterator_tag, ConfigTree>
    {
    public:
        using Iterator = boost::property_tree::ptree::const_assoc_iterator;

        explicit SubtreeIterator(Iterator it, std::string const& root,
                                 ConfigTree const& parent)
            : _it(it), _tagname(root), _parent(parent)
        {}

        SubtreeIterator& operator++() {
            ++_it;
            _has_incremented = true;
            return *this;
        }

        ConfigTree operator*() {
            // if this iterator has been incremented since the last dereference,
            // tell the _parent instance that a subtree now has been parsed.
            if (_has_incremented) {
                _has_incremented = false;
                _parent.markVisited(_tagname, Attr::TAG, false);
            }
            return ConfigTree(_it->second, _parent, _tagname);
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

    protected:
        std::string const _tagname;
        ConfigTree const& _parent;
    };

    /*! A wrapper around a Boost Iterator for iterating over ranges of parameters.
     *
     * The methods of this class tell the associated (parent) \c ConfigTree object when
     * a setting has been parsed.
     */
    class ParameterIterator : public SubtreeIterator
    {
    public:
#if defined(_MSC_VER) && _MSC_VER < 1900
        // 1900 == MSCV 14.0 == Visual Studio 2015
        // according to this post: http://stackoverflow.com/a/70630
        // This table: http://en.cppreference.com/w/cpp/compiler_support
        // says that since MSVC 14.0 inheriting of constructors is supported.
        //! Inherit the constructor
        explicit ParameterIterator(Iterator it, std::string const& root,
                                   ConfigTree const& parent)
            : SubtreeIterator(it, root, parent)
        {}
#else
        //! Inherit the constructor
        using SubtreeIterator::SubtreeIterator;
#endif

        ConfigTree operator*() {
            auto st = SubtreeIterator::operator*();
            if (st.hasChildren())
                _parent.error("The requested parameter <" + _tagname + ">"
                              " has child elements.");
            return st;
        }
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
                               ConfigTree const& parent)
            : _it(it), _tagname(root), _parent(parent)
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
                _parent.markVisited<ValueType>(_tagname, Attr::TAG, false);
            }
            return ConfigTree(_it->second, _parent, _tagname).getValue<ValueType>();
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
        std::string const _tagname;
        ConfigTree const& _parent;
    };

    //! The tree being wrapped by this class.
    using PTree = boost::property_tree::ptree;

    /*! Type of the function objects used as callbacks.
     *
     * Arguments of the callback:
     * \arg \c filename the file being from which this ConfigTree has been read.
     * \arg \c path     the path in the tree where the message was generated.
     * \arg \c message  the message to be printed.
     */
    using Callback = std::function<void(const std::string& filename,
                                        const std::string& path,
                                        const std::string& message)>;

    /*!
     * Creates a new instance wrapping the given Boost Property Tree.
     *
     * \param tree the Boost Property Tree to be wrapped
     * \param filename the file from which the \c tree has been read
     * \param error_cb callback function to be called on error.
     * \param warning_cb callback function to be called on warning.
     *
     * The callback functions must be valid callable functions, i.e. not nullptr's.
     * They are configurable in order to make unit tests of this class easier.
     * They should not be provided in production code!
     *
     * Defaults are strict: By default, both callbacks are set to the same function,
     * i.e., warnings will also result in program abortion!
     */
    explicit ConfigTree(PTree const& tree,
                        std::string const& filename,
                        Callback const& error_cb,
                        Callback const& warning_cb);

    /*! This constructor is deleted in order to prevent the user from passing
     * temporary instances of \c PTree.
     * Doing so would lead to a dangling reference \c _tree and to program crash.
     */
    explicit ConfigTree(PTree&&, std::string const&,
                        Callback const&, Callback const&) = delete;

    //! copying is not compatible with the semantics of this class
    ConfigTree(ConfigTree const&) = delete;

    //! After being moved from, \c other is in an undefined state and must not be
    //! used anymore!
    ConfigTree(ConfigTree && other);

    //! copying is not compatible with the semantics of this class
    ConfigTree& operator=(ConfigTree const&) = delete;

    //! After being moved from, \c other is in an undefined state and must not be
    //! used anymore!
    ConfigTree& operator=(ConfigTree &&);

    /*! \name Methods for directly accessing parameter values
     *
     */
    //!\{

    /*! Get parameter \c param of type \c T from the configuration tree.
     *
     * \return the value looked for.
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    template<typename T> T
    getConfigParameter(std::string const& param) const;

    /*! Get parameter \c param of type \c T from the configuration tree or the \c default_value.
     *
     * This method has a similar behaviour as getConfigParameter(std::string const&) except in case
     * of errors the \c default_value is returned.
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    template<typename T> T
    getConfigParameter(std::string const& param, T const& default_value) const;

    /*! Get parameter \c param of type \c T from the configuration tree if present
     *
     * This method has a similar behaviour as getConfigParameter(std::string const&) except
     * no errors are raised. Rather it can be told from the return value if the
     * parameter could be read.
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    template<typename T> boost::optional<T>
    getConfigParameterOptional(std::string const& param) const;

    /*! Fetches all parameters with name \c param from the current level of the tree.
     *
     * The return value is suitable to be used with range-base for-loops.
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    template<typename T> Range<ValueIterator<T> >
    getConfigParameterList(std::string const& param) const;

    //!\}

    /*! \name Methods for accessing parameters that have attributes
     *
     * The <tt>getConfigParameter...()</tt> methods in this group---note: they do not have template
     * parameters---check that the queried parameters do not have any children (apart from XML
     * attributes); if they do, error() is called.
     *
     * The support for parameters with attributes is limited in the sense that it is not
     * possible to peek/check them. However, such functionality can easily be added on demand.
     */
    //!\{

    /*! Get parameter \c param from the configuration tree.
     *
     * \return the subtree representing the requested parameter
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    ConfigTree
    getConfigParameter(std::string const& param) const;

    /*! Get parameter \c param from the configuration tree if present.
     *
     * \return the subtree representing the requested parameter
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    boost::optional<ConfigTree>
    getConfigParameterOptional(std::string const& param) const;

    /*! Fetches all parameters with name \c param from the current level of the tree.
     *
     * The return value is suitable to be used with range-base for-loops.
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    Range<ParameterIterator>
    getConfigParameterList(std::string const& param) const;

    /*! Get the plain data contained in the current level of the tree.
     *
     * \return the data converted to the type \c T
     *
     * \pre The data must not have been read before.
     */
    template<typename T> T
    getValue() const;

    /*! Get XML attribute \c attr of type \c T for the current parameter.
     *
     * \return the requested attribute
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    template<typename T> T
    getConfigAttribute(std::string const& attr) const;

    /*! Get XML attribute \c attr of type \c T for the current parameter if present.
     *
     * \return the requested attribute
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    template<typename T> boost::optional<T>
    getConfigAttributeOptional(std::string const& attr) const;

    //!\}

    /*! \name Methods for peeking and checking parameters
     *
     * To be used in builder/factory functions: E.g., one can peek a parameter denoting
     * the type of an object to generate in the builder, and check the type parameter in
     * the constructor of the generated object.
     */
    //!\{

    /*! Peek at a parameter \c param of type \c T from the configuration tree.
     *
     * This method is an exception to the single-read rule. It is meant to be used to
     * tell from a ConfigTree instance where to pass that instance on for further processing.
     *
     * But in order that the requested parameter counts as "completely parsed", it has to be
     * read through some other method, too.
     *
     * Return value and error behaviour are the same as for getConfigParameter<T>(std::string const&).
     */
    template<typename T> T
    peekConfigParameter(std::string const& param) const;

    /*! Assert that \c param has the given \c value.
     *
     * Convenience method combining getConfigParameter(std::string const&) with a check.
     */
    template<typename T> void
    checkConfigParameter(std::string const& param, T const& value) const;

    //! Make checkConfigParameter() work for string literals.
    template<typename Ch> void
    checkConfigParameter(std::string const& param, Ch const* value) const;

    //!\}

    /*! \name Methods for accessing subtrees
     *
     */
    //!\{

    /*! Get the subtree rooted at \c root
     *
     * If \c root is not found error() is called.
     *
     * \pre \c root must not have been read before from this ConfigTree.
     */
    ConfigTree
    getConfigSubtree(std::string const& root) const;

    /*! Get the subtree rooted at \c root if present
     *
     * \pre \c root must not have been read before from this ConfigTree.
     */
    boost::optional<ConfigTree>
    getConfigSubtreeOptional(std::string const& root) const;

    /*! Get all subtrees that have a root \c root from the current level of the tree.
     *
     * The return value is suitable to be used with range-base for-loops.
     *
     * \pre \c root must not have been read before from this ConfigTree.
     */
    Range<SubtreeIterator>
    getConfigSubtreeList(std::string const& root) const;

    //!\}

    /*! \name Methods for ignoring parameters
     *
     */
    //!\{

    /*! Tell this instance to ignore parameter \c param.
     *
     * This method is used to avoid warning messages.
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    void ignoreConfigParameter(std::string const& param) const;

    /*! Tell this instance to ignore all parameters \c param on the current level of the tree.
     *
     * This method is used to avoid warning messages.
     *
     * \pre \c param must not have been read before from this ConfigTree.
     */
    void ignoreConfigParameterAll(std::string const& param) const;

    /*! Tell this instance to ignore the XML attribute \c attr.
     *
     * This method is used to avoid warning messages.
     *
     * \pre \c attr must not have been read before from this ConfigTree.
     */
    void ignoreConfigAttribute(std::string const& attr) const;

    //!\}

    //! The destructor performs the check if all nodes at the current level of the tree
    //! have been read.
    //! Errors raised by the check are swallowed. Use assertNoSwallowedErrors()
    //! manually to check for those.
    ~ConfigTree();

    //! Default error callback function
    //! Will throw std::runtime_error
    static void onerror(std::string const& filename, std::string const& path,
                        std::string const& message);

    //! Default warning callback function
    //! Will print a warning message
    static void onwarning(std::string const& filename, std::string const& path,
                          std::string const& message);

    //! Asserts that there have not been any errors reported in the destructor.
    static void assertNoSwallowedErrors();

private:
    //! Default implementation of reading a value of type T.
    template<typename T> boost::optional<T>
    getConfigParameterOptionalImpl(std::string const& param, T*) const;

    //! Implementation of reading a vector of values of type T.
    template<typename T> boost::optional<std::vector<T>>
    getConfigParameterOptionalImpl(std::string const& param, std::vector<T>*) const;

    struct CountType
    {
        int count;
        std::type_index type;
    };

    //! Used to indicate if dealing with XML tags or XML attributes
    enum class Attr : bool
    {
        TAG = false, ATTR = true
    };

    //! Used for wrapping a subtree
    explicit ConfigTree(PTree const& tree, ConfigTree const& parent, std::string const& root);

    /*! Called if an error occurs. Will call the error callback.
     *
     * This method only acts as a helper method and throws std::runtime_error.
     */
#if defined(_MSC_VER) && _MSC_VER < 1900
    // 1900 == MSCV 14.0 == Visual Studio 2015
    // according to this post: http://stackoverflow.com/a/70630
    // This table: http://en.cppreference.com/w/cpp/compiler_support
    // says that since MSVC 14.0 attributes are supported.
    __declspec(noreturn)
#else
    [[noreturn]]
#endif
    void error(std::string const& message) const;

    //! Called for printing warning messages. Will call the warning callback.
    //! This method only acts as a helper method.
    void warning(std::string const& message) const;

    //! Checks if \c key complies with the rules [a-z0-9_].
    void checkKeyname(std::string const& key) const;

    //! Used to generate the path of a subtree.
    std::string joinPaths(std::string const& p1, std::string const& p2) const;

    //! Asserts that the \c key has not been read yet.
    void checkUnique(std::string const& key) const;

    //! Asserts that the attribute \c attr has not been read yet.
    void checkUniqueAttr(std::string const& attr) const;

    /*! Keeps track of the key \c key and its value type \c T.
     *
     * This method asserts that a key is read always with the same type.
     *
     * \c param peek_only if true, do not change the read-count of the given key.
     */
    template<typename T>
    CountType& markVisited(std::string const& key, Attr const is_attr,
                           bool peek_only) const;

    /*! Keeps track of the key \c key and its value type ConfigTree.
     *
     * This method asserts that a key is read always with the same type.
     *
     * \c param peek_only if true, do not change the read-count of the given key.
     */
    CountType& markVisited(std::string const& key, Attr const is_attr,
                           bool const peek_only) const;

    //! Used in the destructor to compute the difference between number of reads of a parameter
    //! and the number of times it exists in the ConfigTree
    void markVisitedDecrement(Attr const is_attr, std::string const& key) const;

    //! Checks if this tree has any children.
    bool hasChildren() const;

    /*! Checks if the top level of this tree has been read entirely (and not too often).
     *
     * \post This method also invalidates the instance, i.e., afterwards it must not
     *       be used anymore!
     */
    void checkAndInvalidate();

    //! returns a short string at suitable for error/warning messages
    static std::string shortString(std::string const& s);

    //! The wrapped tree.
    boost::property_tree::ptree const* _tree;

    //! A path printed in error/warning messages.
    std::string _path;

    //! The path of the file from which this tree has been read.
    std::string _filename;

    //! A pair (is attribute, tag/attribute name).
    using KeyType = std::pair<Attr, std::string>;

    //! A map KeyType -> (count, type) keeping track which parameters have been read
    //! how often and which datatype they have.
    //!
    //! This member will be written to when reading from the config tree.
    //! Therefore it has to be mutable in order to be able to read from
    //! constant instances, e.g., those passed as constant references to
    //! temporaries.
    mutable std::map<KeyType, CountType> _visited_params;

    //! Indicates if the plain data contained in this tree has already been read.
    mutable bool _have_read_data = false;

    Callback _onerror;   //!< Custom error callback.
    Callback _onwarning; //!< Custom warning callback.

    //! Character separating two path components.
    static const char pathseparator;

    //! Set of allowed characters as the first letter of a key name.
    static const std::string key_chars_start;

    //! Set of allowed characters in a key name.
    static const std::string key_chars;

    friend void checkAndInvalidate(ConfigTree* const conf);
    friend void checkAndInvalidate(ConfigTree& conf);
    friend void checkAndInvalidate(std::unique_ptr<ConfigTree> const& conf);
};

}

#include "ConfigTree-impl.h"
