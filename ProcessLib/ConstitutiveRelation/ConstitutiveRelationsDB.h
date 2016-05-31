/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_CONSTITUTIVERELATION_CONSTITUTIVERELATIONSDB_H
#define PROCESSLIB_CONSTITUTIVERELATION_CONSTITUTIVERELATIONSDB_H

#include <iostream>
#include <map>
#include <string>

#include "ConstitutiveRelationBuilder.h"

namespace ProcessLib {
namespace ConstitutiveRelation {

namespace detail
{

template<typename TRef, typename... Ts>
struct AllTypesSameAs;

template<typename TRef, typename T1, typename... Ts>
struct AllTypesSameAs<TRef, T1, Ts...>
{
    static constexpr bool value = std::is_same<TRef, T1>::value
        && AllTypesSameAs<TRef, Ts...>::value;
};

template<typename TRef>
struct AllTypesSameAs<TRef>
{
    static constexpr bool value = true;
};

} // namespace detail


struct AddConstitutiveRelationResult
{
    explicit AddConstitutiveRelationResult(std::string const& name_)
        : name(name_)
    {
#ifndef NDEBUG
        std::cout << "-- AddConstitutiveRelationResult name: " << name_ << "\n";
#endif
    }

    std::string const name;
};

extern bool db_dummy;

class ConstitutiveRelationsDB final
{
public:
    static ConstitutiveRelationsDB& singleton();

    template<typename ConstitutiveRelationReturnType,
             typename... ConstitutiveRelationArguments>
    // AddConstitutiveRelationResult
    bool
    add(std::string const& name,
        std::unique_ptr<ConstitutiveRelationBuilder<
            ConstitutiveRelationReturnType,
            ConstitutiveRelationArguments...>>&&
        constitutive_relation_builder,
        bool const* dummy);

    // can be template-meta-programmed
    using Fct1 = double(*)(double);
    Fct1 get1(std::string const& name);
    using Fct2 = double(*)(double, double);
    Fct2 get2(std::string const& name);

private:
    ConstitutiveRelationsDB() = default;

    std::map<std::string,
             std::unique_ptr<ConstitutiveRelationBuilderBase>>
    _fcts;
};


template<typename ConstitutiveRelationReturnType,
         typename... ConstitutiveRelationArguments>
// AddConstitutiveRelationResult
bool
ConstitutiveRelationsDB::add(
    std::string const& name,
    std::unique_ptr<ConstitutiveRelationBuilder<
        ConstitutiveRelationReturnType,
        ConstitutiveRelationArguments...>>&&
    constitutive_relation_builder,
    bool const* dummy)
{
    static_assert(
            detail::AllTypesSameAs<double, ConstitutiveRelationArguments...>::value,
            "Some of the arguments of the passed function are not of the type double.");

    auto res = _fcts.insert(std::make_pair(
                name,
                std::move(constitutive_relation_builder)));

#ifndef NDEBUG
    std::cout << "== added (true/false: " << res.second << ") fct with name `" << name << "' and "
        << sizeof...(ConstitutiveRelationArguments) << " argument(s).\n";
    std::cout << "== dummy& is " << dummy << std::endl;
#endif

    return res.second;
    // return AddConstitutiveRelationResult(name);
}

} // namespace ConstitutiveRelation
} // namespace ProcessLib

#endif // PROCESSLIB_CONSTITUTIVERELATION_CONSTITUTIVERELATIONSDB_H
