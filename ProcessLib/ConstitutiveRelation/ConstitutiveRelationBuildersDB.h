/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_CONSTITUTIVERELATION_CONSTITUTIVERELATIONBUILDERSDB_H
#define PROCESSLIB_CONSTITUTIVERELATION_CONSTITUTIVERELATIONBUILDERSDB_H

#include <iostream>
#include <map>
#include <string>

#include <logog/include/logog.hpp>

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


class ConstitutiveRelationBuildersDB final
{
public:
    static void initialize();

    template<typename ConstitutiveRelationReturnType,
             typename... ConstitutiveRelationArguments>
    static void
    add(std::string const& name,
        std::unique_ptr<ConstitutiveRelationBuilder<
            ConstitutiveRelationReturnType,
            ConstitutiveRelationArguments...>>&&
        constitutive_relation_builder);

    template<typename ConstitutiveRelationReturnType,
             typename... ConstitutiveRelationArguments>
    static
    ConstitutiveRelationBuilder<
        ConstitutiveRelationReturnType,
        ConstitutiveRelationArguments...> const&
    get(std::string const& name);

    template <typename ConstitutiveRelationReturnType,
              typename... ConstitutiveRelationArguments>
    static std::unique_ptr<ConstitutiveRelation<
        ConstitutiveRelationReturnType, ConstitutiveRelationArguments...>>
    invokeBuilder(BaseLib::ConfigTree const& config,
                  std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const&
                      parameters);

private:
    ConstitutiveRelationBuildersDB() = default;

    static ConstitutiveRelationBuildersDB& singleton();

    std::map<std::string,
             std::unique_ptr<ConstitutiveRelationBuilderBase>>
    _fcts;
};


template<typename ConstitutiveRelationReturnType,
         typename... ConstitutiveRelationArguments>
void
ConstitutiveRelationBuildersDB::add(
    std::string const& name,
    std::unique_ptr<ConstitutiveRelationBuilder<
        ConstitutiveRelationReturnType,
        ConstitutiveRelationArguments...>>&&
    constitutive_relation_builder)
{
    static_assert(
        detail::AllTypesSameAs<double, ConstitutiveRelationArguments...>::value,
        "Some of the arguments of the passed function are not of the type double.");

    auto res = singleton()._fcts.insert(std::make_pair(
                name, std::move(constitutive_relation_builder)));

    if (!res.second) {
        ERR("A function with the name `%s' already exists in the constitutive"
            " relations database.", name.c_str());
        std::abort();
    }

    DBUG("CtiveRelDB: Added a function with name `%s' and %i argument(s).",
         name.c_str(), sizeof...(ConstitutiveRelationArguments));
}


template<typename ConstitutiveRelationReturnType,
         typename... ConstitutiveRelationArguments>
ConstitutiveRelationBuilder<
    ConstitutiveRelationReturnType,
    ConstitutiveRelationArguments...> const&
ConstitutiveRelationBuildersDB::get(
    std::string const& name)
{
    static_assert(
        detail::AllTypesSameAs<double, ConstitutiveRelationArguments...>::value,
        "Some of the arguments of the passed function are not of the type double.");

    auto const& fcts = singleton()._fcts;
    auto it = fcts.find(name);

    if (it == fcts.end()) {
        ERR("A function with the name `%s' has not been found in the constitutive"
            " relations database.", name.c_str());
        std::abort();
    }

    auto const* builder_base = it->second.get();
    auto const* builder_concrete = dynamic_cast<
        ConstitutiveRelationBuilder<
            ConstitutiveRelationReturnType,
                ConstitutiveRelationArguments...> const*>(builder_base);

    if (builder_concrete == nullptr) {
        ERR("The stored and the requested type of the constitutive relation"
            " `%s' differ.", name.c_str());
        std::abort();
    }

    return *builder_concrete;
}

template <typename ConstitutiveRelationReturnType,
          typename... ConstitutiveRelationArguments>
std::unique_ptr<ConstitutiveRelation<ConstitutiveRelationReturnType,
                                     ConstitutiveRelationArguments...>>
ConstitutiveRelationBuildersDB::invokeBuilder(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    auto const& builder =
        get<ConstitutiveRelationReturnType, ConstitutiveRelationArguments...>(
            config.peekConfParam<std::string>("type"));
    return builder.createConstitutiveRelation(config, parameters);
}

} // namespace ConstitutiveRelation
} // namespace ProcessLib

#endif // PROCESSLIB_CONSTITUTIVERELATION_CONSTITUTIVERELATIONBUILDERSDB_H
