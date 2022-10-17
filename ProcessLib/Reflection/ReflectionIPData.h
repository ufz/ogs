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

#include <boost/mp11.hpp>

#include "MathLib/KelvinVector.h"
#include "ReflectionData.h"

namespace ProcessLib::Reflection
{
namespace detail
{
// Used in metaprogramming to check if the type T has a member "reflect".
// Implementation based on https://stackoverflow.com/a/257382
template <typename T>
class HasReflect
{
    template <typename C>
    static char test(decltype(&C::reflect));

    template <typename C>
    static double test(...);

public:
    static constexpr bool value = sizeof(test<T>(nullptr)) == sizeof(char);
};

template <typename T>
struct NumberOfComponents;

template <>
struct NumberOfComponents<double> : std::integral_constant<unsigned, 1>
{
};

template <int N>
struct NumberOfComponents<Eigen::Matrix<double, N, 1, Eigen::ColMajor, N, 1>>
    : std::integral_constant<unsigned, N>
{
};

/** A function object taking a local assembler as its argument and returning a
 * <tt>std::vector\<double\></tt> of some specific "flattened" integration point
 * (IP) data.
 *
 * \tparam Dim the space dimension
 * \tparam Accessor_IPDataVecInLocAsm see below
 * \tparam Accessor_CurrentLevelFromIPDataVecElement see below
 *
 * In OGS IP data is usually stored in the local assembler in the following way:
 *
 * \code
 * struct LocAsm
 * {
 *   std::vector<IPData1> ip_data1;
 *   std::vector<IPData2> ip_data2;
 * };
 * \endcode
 *
 * The types \c IPData1 and \c IPData2 might directly contain the IP data or
 * might have \c struct members who contain the IP data, e.g.:
 *
 * \code
 * struct IPData1
 * {
 *   double scalar1;  // IPData1 directly holds IP data
 *   Eigen::Vector<double, 3> vector1;
 * };
 *
 * struct IPData_Level2
 * {
 *   double scalar2;
 *   Eigen::Vector<double, 6> kelvin2;
 * };
 *
 * struct IPData2
 * {
 *   IPData_Level2 level2;  // IPData2 does not directly hold IP data
 * };
 * \endcode
 *
 * \c Accessor_IPDataVecInLocAsm is a function object with signature
 * <tt>LocAsm const\& -\> std::vector\<IPData\> const\&</tt>.
 *
 * \c Accessor_CurrentLevelFromIPDataVecElement is a function object with
 * signature <tt>IPData const\& -\> double (or Eigen::Vector)</tt>,
 * where \c IPData is the "top level" struct contained in the
 * <tt>std::vector\<IPData\></tt>.
 *
 * I.e. the first accessor takes us from the local assembler to the IP data
 * vector and the second accessor takes us from an IP data vector element to the
 * final IP data of type \c double or <tt>Eigen::Vector</tt>.
 *
 * \note This function object transforms Kelvin vector typed IP data to a
 * ParaView compatible symmetric tensor representation.
 */
template <int Dim, typename Accessor_IPDataVecInLocAsm,
          typename Accessor_CurrentLevelFromIPDataVecElement>
struct GetFlattenedIPDataFromLocAsm
{
    static_assert(!std::is_reference_v<Accessor_IPDataVecInLocAsm>);
    static_assert(
        !std::is_reference_v<Accessor_CurrentLevelFromIPDataVecElement>);

    Accessor_IPDataVecInLocAsm accessor_ip_data_vec_in_loc_asm;
    Accessor_CurrentLevelFromIPDataVecElement
        accessor_current_level_from_ip_data_vec_element;

    template <typename LocAsm>
    std::vector<double> operator()(LocAsm const& loc_asm) const
    {
        using IPDataVector = std::remove_cvref_t<
            std::invoke_result_t<Accessor_IPDataVecInLocAsm, LocAsm const&>>;
        using IPDataVectorElement = typename IPDataVector::value_type;

        // the concrete IP data, e.g. double or Eigen::Vector
        using ConcreteIPData = std::remove_cvref_t<
            std::invoke_result_t<Accessor_CurrentLevelFromIPDataVecElement,
                                 IPDataVectorElement const&>>;

        constexpr unsigned num_comp = NumberOfComponents<ConcreteIPData>::value;

        if constexpr (num_comp == 1)
        {
            auto const& ip_data_vector =
                accessor_ip_data_vec_in_loc_asm(loc_asm);
            auto const num_ips = ip_data_vector.size();

            std::vector<double> result(num_ips);

            // TODO optimization if nothing needs to be copied?
            for (std::size_t ip = 0; ip < num_ips; ++ip)
            {
                auto const& ip_data_vector_element = ip_data_vector[ip];
                auto const& ip_data =
                    accessor_current_level_from_ip_data_vec_element(
                        ip_data_vector_element);

                result[ip] = ip_data;
            }

            return result;
        }
        else if constexpr (num_comp ==
                           MathLib::KelvinVector::kelvin_vector_dimensions(Dim))
        {
            auto const& ip_data_vector =
                accessor_ip_data_vec_in_loc_asm(loc_asm);
            auto const num_ips = ip_data_vector.size();

            std::vector<double> result(num_comp * num_ips);

            for (std::size_t ip = 0; ip < num_ips; ++ip)
            {
                auto const& ip_data_vector_element = ip_data_vector[ip];
                auto const& ip_data =
                    accessor_current_level_from_ip_data_vec_element(
                        ip_data_vector_element);

                auto const converted =
                    MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                        ip_data);

                for (unsigned comp = 0; comp < num_comp; ++comp)
                {
                    result[ip * num_comp + comp] = converted[comp];
                }
            }

            return result;
        }
        else
        {
            auto const& ip_data_vector =
                accessor_ip_data_vec_in_loc_asm(loc_asm);
            auto const num_ips = ip_data_vector.size();

            std::vector<double> result(num_comp * num_ips);

            for (std::size_t ip = 0; ip < num_ips; ++ip)
            {
                auto const& ip_data_vector_element = ip_data_vector[ip];
                auto const& ip_data =
                    accessor_current_level_from_ip_data_vec_element(
                        ip_data_vector_element);

                for (unsigned comp = 0; comp < num_comp; ++comp)
                {
                    result[ip * num_comp + comp] = ip_data[comp];
                }
            }

            return result;
        }
    }
};

// Convenience function for template argument deduction with
// GetFlattenedIPDataFromLocAsm
template <int Dim, typename Accessor_IPDataVecInLocAsm,
          typename Accessor_CurrentLevelFromIPDataVecElement>
GetFlattenedIPDataFromLocAsm<
    Dim, std::remove_cvref_t<Accessor_IPDataVecInLocAsm>,
    std::remove_cvref_t<Accessor_CurrentLevelFromIPDataVecElement>>
getFlattenedIPDataFromLocAsm(
    Accessor_IPDataVecInLocAsm accessor_ip_data_vec_in_loc_asm,
    Accessor_CurrentLevelFromIPDataVecElement
        accessor_current_level_from_ip_data_vec_element)
{
    return {std::forward<Accessor_IPDataVecInLocAsm>(
                accessor_ip_data_vec_in_loc_asm),
            std::forward<Accessor_CurrentLevelFromIPDataVecElement>(
                accessor_current_level_from_ip_data_vec_element)};
}

// Convenience function for template argument deduction with
// GetFlattenedIPDataFromLocAsm. Overload of the function above with less
// arguments.
template <int Dim, typename Accessor_IPDataVecInLocAsm>
auto getFlattenedIPDataFromLocAsm(
    Accessor_IPDataVecInLocAsm&& accessor_ip_data_vec_in_loc_asm)
{
    return getFlattenedIPDataFromLocAsm<Dim>(
        std::forward<Accessor_IPDataVecInLocAsm>(
            accessor_ip_data_vec_in_loc_asm),
        std::identity{});
}

/// Calls the given \c callback for each flattened IP data accessor obtained
/// recursively from the given \c reflection data.
///
/// The \c callback function must take (i) the name of the IP data (\c
/// std::string), (ii) their number of components (\c unsigned) and a flattened
/// IP data accessor of type
/// ProcessLib::Reflection::detail::GetFlattenedIPDataFromLocAsm as arguments.
///
/// \see ProcessLib::Reflection::detail::GetFlattenedIPDataFromLocAsm for
/// details, also on the template parameters.
template <int Dim, typename Callback, typename ReflectionDataTuple,
          typename Accessor_IPDataVecInLocAsm,
          typename Accessor_CurrentLevelFromIPDataVecElement>
void forEachReflectedFlattenedIPDataAccessor(
    Callback const& callback, ReflectionDataTuple const& reflection_data,
    Accessor_IPDataVecInLocAsm const& accessor_ip_data_vec_in_loc_asm,
    Accessor_CurrentLevelFromIPDataVecElement const&
        accessor_current_level_from_ip_data_vec_element)
{
    boost::mp11::tuple_for_each(
        reflection_data,
        [&accessor_ip_data_vec_in_loc_asm,
         &accessor_current_level_from_ip_data_vec_element,
         &callback]<typename Class, typename Member>(
            ReflectionData<Class, Member> const& refl_data)
        {
            auto accessor_member_from_ip_data_vec_element =
                [level = refl_data.field,
                 accessor_current_level_from_ip_data_vec_element](
                    auto const& ip_data_vec_element) -> Member const&
            {
                return accessor_current_level_from_ip_data_vec_element(
                           ip_data_vec_element).*
                       level;
            };

            if constexpr (HasReflect<Member>::value)
            {
                forEachReflectedFlattenedIPDataAccessor<Dim>(
                    callback, Member::reflect(),
                    accessor_ip_data_vec_in_loc_asm,
                    accessor_member_from_ip_data_vec_element);
            }
            else
            {
                constexpr unsigned num_comp = NumberOfComponents<Member>::value;

                callback(refl_data.name, num_comp,
                         getFlattenedIPDataFromLocAsm<Dim>(
                             accessor_ip_data_vec_in_loc_asm,
                             accessor_member_from_ip_data_vec_element));
            }
        });
}

// Overload of the function above with less arguments
template <int Dim, typename Callback, typename ReflectionDataTuple,
          typename Accessor_IPDataVecInLocAsm>
void forEachReflectedFlattenedIPDataAccessor(
    Callback const& callback,
    ReflectionDataTuple const& reflection_data,
    Accessor_IPDataVecInLocAsm const& accessor_ip_data_vec_in_loc_asm)
{
    forEachReflectedFlattenedIPDataAccessor<Dim>(
        callback, reflection_data, accessor_ip_data_vec_in_loc_asm,
        std::identity{});
}

}  // namespace detail

/// Calls the passed \c callback for each IP data accessor for the given
/// \c LocAsmIF class.
///
/// The IP data accessors are obtained via reflection (i.e., via the static
/// \c reflect() method.
///
/// The IP data accessors provide IP data as a flat \c std::vector<double>.
///
/// The \c callback must accept name, number of components and a function object
/// with signature <tt>LocAsmIF const\& -\> std::vector\<double\></tt> as its
/// arguments.
template <int Dim, typename LocAsmIF, typename Callback, typename ReflData>
void forEachReflectedFlattenedIPDataAccessor(ReflData const& reflection_data,
                                             Callback const& callback)
{
    boost::mp11::tuple_for_each(
        reflection_data,
        [&callback]<typename Class, typename Member>(
            ReflectionData<Class, std::vector<Member>> const& refl_data)
        {
            auto accessor_ip_data_vec_in_loc_asm = [ip_data_vector =
                                                        refl_data.field](
                LocAsmIF const& loc_asm) -> auto const&
            {
                return loc_asm.*ip_data_vector;
            };

            if constexpr (detail::HasReflect<Member>::value)
            {
                detail::forEachReflectedFlattenedIPDataAccessor<Dim>(
                    callback, Member::reflect(),
                    accessor_ip_data_vec_in_loc_asm);
            }
            else
            {
                constexpr unsigned num_comp =
                    detail::NumberOfComponents<Member>::value;

                callback(refl_data.name, num_comp,
                         detail::getFlattenedIPDataFromLocAsm<Dim>(
                             accessor_ip_data_vec_in_loc_asm));
            }
        });
}
}  // namespace ProcessLib::Reflection
