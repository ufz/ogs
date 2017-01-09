/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>
#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace detail
{
//! Checks if all types \c Ts are the same as \c TRef.
template<typename TRef, typename... Ts>
struct AllTypesSameAs;

template<typename TRef, typename T1, typename... Ts>
struct AllTypesSameAs<TRef, T1, Ts...>
{
    static const bool value = std::is_same<TRef, T1>::value
        && AllTypesSameAs<TRef, Ts...>::value;
};

template<typename TRef>
struct AllTypesSameAs<TRef>
{
    static const bool value = true;
};

template<typename TRef, typename T1, typename... Ts>
const bool AllTypesSameAs<TRef, T1, Ts...>::value;

template<typename TRef>
const bool AllTypesSameAs<TRef>::value;

} // namespace detail

namespace NumLib
{
//! Stores a function object along with a name for it and information about its
//! arguments.
class NamedFunction final
{
public:
    /*! Constructs a new named function.
     *
     * \param name the function's name
     * \param argument_names names  of arguments of the function
     * \param function the actual function object
     */
    template <typename ReturnType, typename... Arguments>
    NamedFunction(std::string const& name,
                  std::vector<std::string>&& argument_names,
                  std::function<ReturnType(Arguments...)>&& function);

    NamedFunction(NamedFunction&& other);
    NamedFunction(NamedFunction const&);

    ~NamedFunction();

    //! Returns the function's name.
    std::string const& getName() const { return _name; }
    //! Returns the names of the function's arguments.
    std::vector<std::string> const& getArgumentNames() const
    {
        return _argument_names;
    }

    //! Call the function with the supplied arguments.
    double call(std::vector<double> const& arguments) const;

    //! Maximum number of function arguments supported by NamedFunction.
    static const int MAX_FUNCTION_ARGS = 32;

private:
    //! The function's name.
    std::string _name;

    //! Information about the function's arguments.
    std::vector<std::string> _argument_names;

    //! The function handle.
    void* _function;
};

template <typename ReturnType, typename... Arguments>
NamedFunction::NamedFunction(std::string const& name,
                             std::vector<std::string>&& argument_names,
                             std::function<ReturnType(Arguments...)>&& function)
    : _name(name),
      _argument_names(std::move(argument_names)),
      _function(
          new std::function<ReturnType(Arguments...)>(std::move(function)))
{
    static_assert(
        ::detail::AllTypesSameAs<double, ReturnType, Arguments...>::value,
        "Some of the arguments or the return type of the passed function are "
        "not of the type double.");
    static_assert(sizeof...(Arguments) <= MAX_FUNCTION_ARGS,
                  "The function you passed has too many arguments.");

    assert(sizeof...(Arguments) == _argument_names.size());
}

}  // namespace NumLib
