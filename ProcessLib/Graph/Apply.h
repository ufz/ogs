/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Get.h"

namespace ProcessLib::Graph
{
namespace detail
{
template <typename Function>
struct GetFunctionArgumentTypesPlain  // plain, i.e., without cvref
    /** \cond */
    : GetFunctionArgumentTypesPlain<decltype(&Function::operator())>
/** \endcond */
{
};

// member function
template <typename Result, typename Class, typename... Args>
struct GetFunctionArgumentTypesPlain<Result (Class::*)(Args...)>
{
    using type = boost::mp11::mp_list<std::remove_cvref_t<Args>...>;
};

// const member function
template <typename Result, typename Class, typename... Args>
struct GetFunctionArgumentTypesPlain<Result (Class::*)(Args...) const>
{
    using type = boost::mp11::mp_list<std::remove_cvref_t<Args>...>;
};

// standalone function
template <typename Result, typename... Args>
struct GetFunctionArgumentTypesPlain<Result (*)(Args...)>
{
    using type = boost::mp11::mp_list<std::remove_cvref_t<Args>...>;
};

template <typename Function>
struct GetFunctionReturnType
    /** \cond */
    : GetFunctionReturnType<decltype(&Function::operator())>
/** \endcond */
{
};

// member function
template <typename Result, typename Class, typename... Args>
struct GetFunctionReturnType<Result (Class::*)(Args...)>
{
    using type = Result;
};

// const member function
template <typename Result, typename Class, typename... Args>
struct GetFunctionReturnType<Result (Class::*)(Args...) const>
{
    using type = Result;
};

// standalone function
template <typename Result, typename... Args>
struct GetFunctionReturnType<Result (*)(Args...)>
{
    using type = Result;
};

// Invokes the passed function object with arguments taken ("unpacked") from the
// passed tuples ts. The arguments are picked/passed according to their data
// type.
//
// overload for function objects
template <typename Object,
          typename... Tuples,
          typename... MemberFunctionArgumentTypesPlain>
auto unpackAndInvoke(boost::mp11::mp_list<MemberFunctionArgumentTypesPlain...>,
                     Object&& o,
                     Tuples&... ts) ->
    // std::decay_t "decays" lambdas
    typename GetFunctionReturnType<
        decltype(&std::decay_t<Object>::operator())>::type
{
    return o(
        ProcessLib::Graph::get<MemberFunctionArgumentTypesPlain>(ts...)...);
}

// member function
template <typename Result,
          typename Object,
          typename... Args,
          typename... Tuples,
          typename... MemberFunctionArgumentTypesPlain>
Result unpackAndInvoke(
    boost::mp11::mp_list<MemberFunctionArgumentTypesPlain...>,
    Result (Object::*m)(Args...),
    Object& o,
    Tuples&... ts)
{
    return (o.*m)(
        ProcessLib::Graph::get<MemberFunctionArgumentTypesPlain>(ts...)...);
}

// member function & temporary object
template <typename Result,
          typename Object,
          typename... Args,
          typename... Tuples,
          typename... MemberFunctionArgumentTypesPlain>
Result unpackAndInvoke(
    boost::mp11::mp_list<MemberFunctionArgumentTypesPlain...>,
    Result (Object::*m)(Args...),
    Object&& o,
    Tuples&... ts)
{
    return (o.*m)(
        ProcessLib::Graph::get<MemberFunctionArgumentTypesPlain>(ts...)...);
}

// const member function
template <typename Result,
          typename Object,
          typename... Args,
          typename... Tuples,
          typename... MemberFunctionArgumentTypesPlain>
Result unpackAndInvoke(
    boost::mp11::mp_list<MemberFunctionArgumentTypesPlain...>,
    Result (Object::*m)(Args...) const,
    Object const& o,
    Tuples&... ts)
{
    return (o.*m)(
        ProcessLib::Graph::get<MemberFunctionArgumentTypesPlain>(ts...)...);
}

// standalone function
template <typename Result,
          typename... Args,
          typename... Tuples,
          typename... FunctionArgumentTypesPlain>
Result unpackAndInvoke(boost::mp11::mp_list<FunctionArgumentTypesPlain...>,
                       Result (*fct)(Args...),
                       Tuples&... ts)
{
    return fct(ProcessLib::Graph::get<FunctionArgumentTypesPlain>(ts...)...);
}

template <typename Function, typename... PassedArgsTuples>
struct GetFlattenedTupleOfPassedArgs
    /** \cond */
    : GetFlattenedTupleOfPassedArgs<decltype(&Function::operator()),
                                    Function,
                                    PassedArgsTuples...>
/** \endcond */
{
};

// member function
template <typename Result,
          typename Class,
          typename... DeclaredArgs,
          typename ClassAgain,
          typename... PassedArgsTuples>
struct GetFlattenedTupleOfPassedArgs<Result (Class::*)(DeclaredArgs...),
                                     // first "argument" is the object to which
                                     // the member function belongs
                                     ClassAgain,
                                     PassedArgsTuples...>
{
    using type = typename detail::GetFlattenedTupleTypes<
        std::remove_cvref_t<PassedArgsTuples>...>::type;
};

// const member function
template <typename Result,
          typename Class,
          typename... DeclaredArgs,
          typename... ClassAndPassedArgsTuples>
struct GetFlattenedTupleOfPassedArgs<Result (Class::*)(DeclaredArgs...) const,
                                     ClassAndPassedArgsTuples...>
    // redirect to non-const member function implementation
    : GetFlattenedTupleOfPassedArgs<Result (Class::*)(DeclaredArgs...),
                                    ClassAndPassedArgsTuples...>
{
};

// standalone function
template <typename Result, typename... DeclaredArgs, typename... Tuples>
struct GetFlattenedTupleOfPassedArgs<Result (*)(DeclaredArgs...), Tuples...>
{
    using type = typename detail::GetFlattenedTupleTypes<
        std::remove_cvref_t<Tuples>...>::type;
};

/**
 * Invokes the passed function \c f.
 *
 * \c f can be a function object, a free function, or a member function pointer.
 * In the latter case, the first entry of \c args must be the object for which
 * \c f shall be invoked. (The remainder of) \c args are an arbitrary number of
 * tuples from which arguments to \c f are taken based on their types.
 *
 * \see ProcessLib::Graph::apply()
 */
template <typename Function, typename... Args>
auto applyImpl(Function&& f, Args&&... args) ->
    typename detail::GetFunctionReturnType<std::decay_t<Function>>::type
{
    using namespace boost::mp11;

    using FunctionPlain = std::decay_t<Function>;
    using FlattenedTuple =
        typename GetFlattenedTupleOfPassedArgs<FunctionPlain, Args...>::type;
    using FlattenedTupleOfPlainTypes =
        mp_transform<std::remove_cvref_t, FlattenedTuple>;

    static_assert(
        boost::mp11::mp_is_set<FlattenedTupleOfPlainTypes>::value,
        "The types of all elements of all passed tuples must be unique.");

    using FunctionArgumentTypesPlain =
        typename detail::GetFunctionArgumentTypesPlain<FunctionPlain>::type;

    static_assert(
        boost::mp11::mp_is_set<FunctionArgumentTypesPlain>::value,
        "The argument types of the function to be called must be unique.");

    return unpackAndInvoke(FunctionArgumentTypesPlain{},
                           std::forward<Function>(f),
                           std::forward<Args>(args)...);
}
}  // namespace detail

/**
 * Invokes the passed function (object) \c f with arguments taken from the
 * passed tuples.
 *
 * The passed arguments are determined from their types. Therefore, both the
 * argument types of the function and the member types of all passed tuples must
 * be unique.
 */
template <typename Function, typename... Tuples>
auto apply(Function& f, Tuples&... ts) ->
    typename detail::GetFunctionReturnType<std::decay_t<Function>>::type
{
    static_assert(boost::mp11::mp_similar<std::tuple<>,
                                          std::remove_cv_t<Tuples>...>::value,
                  "In the argument ts... there must be only std::tuple's.");

    return detail::applyImpl(f, ts...);
}

/**
 * Invokes the eval() method of the passed object \c f with arguments taken from
 * the passed tuples.
 *
 * \see apply()
 */
template <typename Function, typename... Tuples>
auto eval(Function& f, Tuples&... ts) ->
    typename detail::GetFunctionReturnType<decltype(&Function::eval)>::type
{
    static_assert(boost::mp11::mp_similar<std::tuple<>,
                                          std::remove_cv_t<Tuples>...>::value,
                  "In the argument ts... there must be only std::tuple's.");

    return detail::applyImpl(&Function::eval, f, ts...);
}
}  // namespace ProcessLib::Graph
