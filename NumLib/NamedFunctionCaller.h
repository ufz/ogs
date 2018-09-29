/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <vector>
#include "NamedFunction.h"

namespace NumLib
{
class SpecificFunctionCaller;

//! Builds expression trees of named functions dynamically at runtime.
class NamedFunctionCaller final
{
public:
    //! Constructs an instance whose unbound arguments have the given names.
    explicit NamedFunctionCaller(
        std::initializer_list<std::string> unbound_argument_names);

    //! Adds the given named function
    void addNamedFunction(NamedFunction&& fct);

    //! Returns all named functions associated with the caller instance.
    std::vector<NamedFunction> const& getNamedFunctions() const
    {
        return _named_functions;
    }

    //! Declares that the argument with name \c sink_arg of the function \c
    //! sink_fct is being computed by the function \c source_fct.
    //!
    //! The functions involved need not already be known to the
    //! NamedFunctionCaller.
    void plug(std::string const& sink_fct, std::string const& sink_arg,
              std::string const& source_fct);

    //! Actually plug all the plugs previously declared.
    //!
    //! \pre All functions involved must have been added.
    void applyPlugs();

    //! Creates a function caller that is able to call the function with the
    //! given name.
    SpecificFunctionCaller getSpecificFunctionCaller(
        std::string const& function_name);

    //! Returns a string representing the expression graph of the given
    //! function.
    //!
    //! \pre applyPlugs() must have been called before.
    std::string getCallExpression(std::string const& function_name) const;

    //! Returns the number of unbound arguments.
    std::size_t getNumberOfUnboundArguments() const;

private:
    //! Calls the function with the given index with the given unbound
    //! arguments.
    double call(std::size_t function_idx,
                std::vector<double> const& unbound_arguments) const;

    //! Maps function names to indices.
    //! Negative indices refer to unbound arguments.
    std::map<std::string, int> _map_name_idx;

    //! Contains all named functions.
    std::vector<NamedFunction> _named_functions;

    //! The expression graph.
    //! Contains for each named function (outer vector) a vector which maps each
    //! function argument to the source function index that computes this
    //! argument.
    std::vector<std::vector<int>> _map_sink_source;

    //! Magic number used to mark function arguments in \c _map_sink_source
    //! whose source functions have not yet been set up.
    const int _uninitialized;

    struct SinkSource
    {
        std::string const sink_fct;
        std::string const sink_arg;
        std::string const source;
    };

    //! Saves plugs declared by plug().
    std::vector<SinkSource> _deferred_plugs;

    friend class SpecificFunctionCaller;
};

//! A function caller that can call one specific function.
//!
//! \todo Use this class to provide some optimizations of the expression
//! evaluation.
class SpecificFunctionCaller final
{
public:
    //! Constructs a new instance.
    SpecificFunctionCaller(std::size_t const function_idx,
                          NamedFunctionCaller const& caller);

    //! Call the function set up with the given unbound arguments.
    double call(std::vector<double> const& unbound_arguments) const;

    //! Returns the number of unbound arguments.
    std::size_t getNumberOfUnboundArguments() const;

private:
    //! Index of the referenced function.
    std::size_t const _function_idx;

    //! The named function caller used for the evaluation.
    NamedFunctionCaller const& _caller;
};

} // namespace NumLib
