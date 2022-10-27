/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Core>
#include <vector>

#include "MathLib/KelvinVector.h"
#include "ReflectionIPData.h"

namespace ProcessLib::Reflection
{
namespace detail
{
template <int dim, typename IPData, typename Accessor>
void setIPData(double const* values,
               std::vector<IPData>& ip_data_vector,
               Accessor const& accessor)
{
    using AccessorResult = std::invoke_result_t<Accessor, IPData&>;
    using AccessorResultStripped = std::remove_cvref_t<AccessorResult>;
    constexpr auto num_comp = NumberOfComponents<AccessorResultStripped>::value;

    auto constexpr kv_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(dim);
    auto const num_int_pts = ip_data_vector.size();

    if constexpr (num_comp == 1)
    {
        for (unsigned ip = 0; ip < num_int_pts; ++ip)
        {
            accessor(ip_data_vector[ip]) = values[ip];
        }
    }
    else
    {
        auto const values_mat =
            Eigen::Map<Eigen::Matrix<double, num_comp, Eigen::Dynamic,
                                     Eigen::ColMajor> const>(values, num_comp,
                                                             num_int_pts);

        for (unsigned ip = 0; ip < num_int_pts; ++ip)
        {
            if constexpr (num_comp == kv_size)
            {
                accessor(ip_data_vector[ip]) =
                    MathLib::KelvinVector::symmetricTensorToKelvinVector(
                        values_mat.col(ip));
            }
            else
            {
                accessor(ip_data_vector[ip]) = values_mat.col(ip);
            }
        }
    }
}

// Sets IP data if the passed name equals the one in refl_data.
// Returns true if IP have been set, false otherwise.
template <int dim, typename IPData, typename Accessor_CurrentLevelFromIPData,
          typename Class, typename Member>
bool setIPDataIfNameMatches(std::string const& name, double const* values,
                            std::vector<IPData>& ip_data_vector,
                            Accessor_CurrentLevelFromIPData const& accessor,
                            ReflectionData<Class, Member> const& refl_data)
{
    auto const field = refl_data.field;

    auto const accessor_field_from_ip_data = [accessor,
                                              field](IPData& ip_data) -> Member&
    { return accessor(ip_data).*field; };

    if constexpr (detail::has_reflect<Member>)
    {
        return reflectSetIPData<dim>(name, values, ip_data_vector,
                                     accessor_field_from_ip_data,
                                     Member::reflect());
    }
    else
    {
        if (refl_data.name != name)
        {
            return false;
        }

        setIPData<dim>(values, ip_data_vector, accessor_field_from_ip_data);

        return true;
    }
}

template <int dim, typename IPData, typename Accessor_CurrentLevelFromIPData,
          typename... Class, typename... Member, std::size_t... Idcs>
bool reflectSetIPData(
    std::string const& name, double const* values,
    std::vector<IPData>& ip_data_vector,
    Accessor_CurrentLevelFromIPData const& accessor,
    std::tuple<ReflectionData<Class, Member>...> const& refl_data,
    std::index_sequence<Idcs...>)
{
    // uses short-circuit evaluation of the fold || ... to stop after the first
    // match
    return ((setIPDataIfNameMatches<dim>(name, values, ip_data_vector, accessor,
                                         std::get<Idcs>(refl_data))) ||
            ...);
}

template <int dim, typename IPData, typename Accessor_CurrentLevelFromIPData,
          typename... Class, typename... Member>
bool reflectSetIPData(
    std::string const& name, double const* values,
    std::vector<IPData>& ip_data_vector,
    Accessor_CurrentLevelFromIPData const& accessor,
    std::tuple<ReflectionData<Class, Member>...> const& refl_data)
{
    return reflectSetIPData<dim>(name, values, ip_data_vector, accessor,
                                 refl_data,
                                 std::make_index_sequence<sizeof...(Class)>{});
}
}  // namespace detail

/**
 * Sets integration point data for the property with the given \c name to the
 * passed \c values.
 *
 * Possible candidate properties are obtained from \c IPData via some sort of
 * reflection.
 *
 * \return The number of integration points.
 */
template <int dim, typename IPData>
std::size_t reflectSetIPData(std::string const& name, double const* values,
                             std::vector<IPData>& ip_data_vector)
{
    detail::reflectSetIPData<dim>(name, values, ip_data_vector, std::identity{},
                                  IPData::reflect());

    return ip_data_vector.size();
}

}  // namespace ProcessLib::Reflection
