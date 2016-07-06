/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_NAMED_FUNCTION
#define NUMLIB_NAMED_FUNCTION

#include <string>
#include <vector>
#include <functional>
#include <memory>
#include <cassert>
#include <type_traits>

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
//! Maximum number of function arguments supported by NamedFunction.
const int MAX_FUNCTION_ARGS = 32;

//! Stores a function object along with a name for it and information about its
//! arguments.
class NamedFunction final
{
public:
    using ArgumentInfo = std::string;

    /*! Constructs a new named function.
     *
     * \param name the function's name
     * \param arguments names  of arguments of the function
     * \param function the actual function object
     */
    template <typename ReturnType, typename... Arguments>
    NamedFunction(
        std::string const& name,
        std::vector<ArgumentInfo>&& arguments,
        std::function<ReturnType(Arguments...)>&& function);

    NamedFunction(NamedFunction&& other);
    NamedFunction(NamedFunction const&);

    ~NamedFunction();

    //! Returns the function's name.
    std::string const& getName() const { return _name; }

    //! Returns information about the function's arguments.
    std::vector<ArgumentInfo> const& getArgumentInfo() const
    {
        return _argument_info;
    }

    //! Call the function with the supplied arguments.
    double call(std::vector<double> const& arguments) const;

private:
    //! The function's name.
    std::string _name;

    //! Information about the function's arguments.
    std::vector<ArgumentInfo> _argument_info;

    //! The function handle.
    void* _function;
};

template <typename ReturnType, typename... Arguments>
NamedFunction::NamedFunction(std::string const& name,
                             std::vector<ArgumentInfo>&& argument_info,
                             std::function<ReturnType(Arguments...)>&& function)
    : _name(name),
      _argument_info(std::move(argument_info)),
      _function(
          new std::function<ReturnType(Arguments...)>(std::move(function)))
{
    static_assert(
        ::detail::AllTypesSameAs<double, ReturnType, Arguments...>::value,
        "Some of the arguments or the return type of the passed function are "
        "not of the type double.");
    static_assert(sizeof...(Arguments) <= MAX_FUNCTION_ARGS,
                  "The function you passed has too many arguments.");

    assert(sizeof...(Arguments) == _argument_info.size());
}

} // namespace NumLib


#endif // NUMLIB_NAMED_FUNCTION
