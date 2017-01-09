/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NamedFunction.h"
#include "BaseLib/TMPUtil.h"

//! Helper struct used in conjunction with BaseLib::IntegerSequence to get a
//! sequence of types, where each type is double.
template <int>
struct Double
{
    using type = double;
};

/*! Calls the given \c function with the given  \c arguments.
 *
 * \tparam Indices sequence of integers used to expand the \c arguments vector
 * to individual arguments.
 */
template <int... Indices>
double call_(void* function, std::vector<double> const& arguments)
{
    assert(arguments.size() == sizeof...(Indices));
    auto* fct = reinterpret_cast<
        std::function<double(typename Double<Indices>::type...)>*>(function);
    auto* args = arguments.data();
    return (*fct)(args[Indices]...);
}

typedef double (*CallerFunction) (void*, std::vector<double> const&);

//! Helps instantiating the call_() function.
template <int... Indices>
CallerFunction generateCallerFromIntegerSequence(
    BaseLib::IntegerSequence<Indices...>)
{
    return call_<Indices...>;
}

//! Instantiates the call_() function for the provided number of arguments.
template <int NumArguments>
CallerFunction generateCaller()
{
    return generateCallerFromIntegerSequence(
        typename BaseLib::GenerateIntegerSequence<NumArguments>::type{});
}

//! Holds instantiations of the call_() function for various numbers of
//! arguments.
static const CallerFunction callers[] = {
    generateCaller<0>(),  generateCaller<1>(),  generateCaller<2>(),
    generateCaller<3>(),  generateCaller<4>(),  generateCaller<5>(),
    generateCaller<6>(),  generateCaller<7>(),  generateCaller<8>(),
    generateCaller<9>(),  generateCaller<10>(), generateCaller<11>(),
    generateCaller<12>(), generateCaller<13>(), generateCaller<14>(),
    generateCaller<15>(), generateCaller<16>(), generateCaller<17>(),
    generateCaller<18>(), generateCaller<19>(), generateCaller<20>(),
    generateCaller<21>(), generateCaller<22>(), generateCaller<23>(),
    generateCaller<24>(), generateCaller<25>(), generateCaller<26>(),
    generateCaller<27>(), generateCaller<28>(), generateCaller<29>(),
    generateCaller<30>(), generateCaller<31>(), generateCaller<32>()};
static_assert(sizeof(callers) / sizeof(CallerFunction) ==
                  NumLib::NamedFunction::MAX_FUNCTION_ARGS + 1,
              "You did not instantiate the right number of callers.");

/*! Deletes the given \c function.
 *
 * \tparam Indices sequence of integers used to cast the \c function to the
 * correct type.
 */
template <int... Indices>
void delete_(void* function)
{
    auto* fct = reinterpret_cast<
        std::function<double(typename Double<Indices>::type...)>*>(function);
    delete fct;
}

typedef void (*DeleterFunction) (void*);

//! Helps instantiating the delete_() function.
template <int... Indices>
DeleterFunction generateDeleterFromIntegerSequence(
    BaseLib::IntegerSequence<Indices...>)
{
    return delete_<Indices...>;
}

//! Instantiates the delete_() function for the provided number of arguments.
template <int NumArguments>
DeleterFunction generateDeleter()
{
    return generateDeleterFromIntegerSequence(
        typename BaseLib::GenerateIntegerSequence<NumArguments>::type{});
}

//! Holds instantiations of the delete_() function for various numbers of
//! arguments.
static const DeleterFunction deleters[] = {
    generateDeleter<0>(),  generateDeleter<1>(),  generateDeleter<2>(),
    generateDeleter<3>(),  generateDeleter<4>(),  generateDeleter<5>(),
    generateDeleter<6>(),  generateDeleter<7>(),  generateDeleter<8>(),
    generateDeleter<9>(),  generateDeleter<10>(), generateDeleter<11>(),
    generateDeleter<12>(), generateDeleter<13>(), generateDeleter<14>(),
    generateDeleter<15>(), generateDeleter<16>(), generateDeleter<17>(),
    generateDeleter<18>(), generateDeleter<19>(), generateDeleter<20>(),
    generateDeleter<21>(), generateDeleter<22>(), generateDeleter<23>(),
    generateDeleter<24>(), generateDeleter<25>(), generateDeleter<26>(),
    generateDeleter<27>(), generateDeleter<28>(), generateDeleter<29>(),
    generateDeleter<30>(), generateDeleter<31>(), generateDeleter<32>()};
static_assert(sizeof(deleters) / sizeof(DeleterFunction) ==
                  NumLib::NamedFunction::MAX_FUNCTION_ARGS + 1,
              "You did not instantiate the right number of deleters.");

/*! Copies the given \c function.
 *
 * \tparam Indices sequence of integers used to cast the \c function to the
 * correct type.
 */
template <int... Indices>
void* copy_(void* function)
{
    auto* fct = reinterpret_cast<
        std::function<double(typename Double<Indices>::type...)>*>(function);
    return new std::function<double(typename Double<Indices>::type...)>(*fct);
}

typedef void* (*CopierFunction) (void*);

//! Helps instantiating the copy_() function.
template <int... Indices>
CopierFunction generateCopierFromIntegerSequence(
    BaseLib::IntegerSequence<Indices...>)
{
    return copy_<Indices...>;
}

//! Instantiates the copy_() function for the provided number of arguments.
template <int NumArguments>
CopierFunction generateCopier()
{
    return generateCopierFromIntegerSequence(
        typename BaseLib::GenerateIntegerSequence<NumArguments>::type{});
}

//! Holds instantiations of the copy_() function for various numbers of
//! arguments.
static const CopierFunction copiers[] = {
    generateCopier<0>(),  generateCopier<1>(),  generateCopier<2>(),
    generateCopier<3>(),  generateCopier<4>(),  generateCopier<5>(),
    generateCopier<6>(),  generateCopier<7>(),  generateCopier<8>(),
    generateCopier<9>(),  generateCopier<10>(), generateCopier<11>(),
    generateCopier<12>(), generateCopier<13>(), generateCopier<14>(),
    generateCopier<15>(), generateCopier<16>(), generateCopier<17>(),
    generateCopier<18>(), generateCopier<19>(), generateCopier<20>(),
    generateCopier<21>(), generateCopier<22>(), generateCopier<23>(),
    generateCopier<24>(), generateCopier<25>(), generateCopier<26>(),
    generateCopier<27>(), generateCopier<28>(), generateCopier<29>(),
    generateCopier<30>(), generateCopier<31>(), generateCopier<32>()};
static_assert(sizeof(copiers) / sizeof(CopierFunction) ==
                  NumLib::NamedFunction::MAX_FUNCTION_ARGS + 1,
              "You did not instantiate the right number of deleters.");

namespace NumLib
{
NamedFunction::NamedFunction(NamedFunction&& other)
    : _name(std::move(other._name)),
      _argument_names(std::move(other._argument_names)),
      _function(other._function)
{
    other._function = nullptr;
}

NamedFunction::NamedFunction(NamedFunction const& other)
    : _name(other._name),
      _argument_names(other._argument_names),
      _function(copiers[_argument_names.size()](other._function))
{
}

NamedFunction::~NamedFunction()
{
    deleters[_argument_names.size()](_function);
}

double NamedFunction::call(const std::vector<double>& arguments) const
{
    assert(arguments.size() == _argument_names.size());
    return callers[_argument_names.size()](_function, arguments);
}

}  // namespace NumLib
