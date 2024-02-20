/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
template <typename T>
concept has_reflect = requires { T::reflect(); };

template <typename... Ts>
auto reflect(std::type_identity<std::tuple<Ts...>>)
{
    using namespace boost::mp11;

    // The types Ts... must be unique. Duplicate types are incompatible with the
    // concept of "reflected" I/O: they would lead to duplicate names for the
    // I/O data.
    static_assert(mp_is_set<mp_list<Ts...>>::value);

    return reflectWithoutName<std::tuple<Ts...>>(
        [](auto& tuple_) -> auto& { return std::get<Ts>(tuple_); }...);
}

template <has_reflect T>
auto reflect(std::type_identity<T>)
{
    return T::reflect();
}

template <typename T>
concept is_reflectable = requires {
    ProcessLib::Reflection::detail::reflect(std::type_identity<T>{});
};

/**
 * Raw data is data that will be read or written, e.g., double values or Eigen
 * vectors.
 *
 * Non-raw data is data for which further reflection will be performed (to find
 * out the members).
 */
template <typename T>
struct is_raw_data : std::false_type
{
};

template <>
struct is_raw_data<double> : std::true_type
{
};

template <int N>
struct is_raw_data<Eigen::Matrix<double, N, 1, Eigen::ColMajor, N, 1>>
    : std::true_type
{
};

template <int N, int M>
struct is_raw_data<Eigen::Matrix<double, N, M, Eigen::RowMajor, N, M>>
    : std::true_type
{
};

template <typename T>
constexpr bool is_raw_data_v = is_raw_data<T>::value;

template <typename T>
struct NumberOfRows;

template <>
struct NumberOfRows<double> : std::integral_constant<unsigned, 1>
{
};

template <int N>
struct NumberOfRows<Eigen::Matrix<double, N, 1, Eigen::ColMajor, N, 1>>
    : std::integral_constant<unsigned, N>
{
};

template <int N, int M>
struct NumberOfRows<Eigen::Matrix<double, N, M, Eigen::RowMajor, N, M>>
    : std::integral_constant<unsigned, N>
{
};

template <typename T>
struct NumberOfColumns;

template <>
struct NumberOfColumns<double> : std::integral_constant<unsigned, 1>
{
};

template <int N>
struct NumberOfColumns<Eigen::Matrix<double, N, 1, Eigen::ColMajor, N, 1>>
    : std::integral_constant<unsigned, 1>
{
};

template <int N, int M>
struct NumberOfColumns<Eigen::Matrix<double, N, M, Eigen::RowMajor, N, M>>
    : std::integral_constant<unsigned, M>
{
};

template <typename T>
struct NumberOfComponents
    : std::integral_constant<unsigned,
                             NumberOfRows<T>::value * NumberOfColumns<T>::value>
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
        static_assert(is_raw_data<ConcreteIPData>::value,
                      "This method only deals with raw data. The given "
                      "ConcreteIPData is not raw data.");

        constexpr unsigned num_rows = NumberOfRows<ConcreteIPData>::value;
        constexpr unsigned num_cols = NumberOfColumns<ConcreteIPData>::value;
        constexpr unsigned num_comp = num_rows * num_cols;
        auto const& ip_data_vector = accessor_ip_data_vec_in_loc_asm(loc_asm);
        auto const num_ips = ip_data_vector.size();

        std::vector<double> result(num_comp * num_ips);

        for (std::size_t ip = 0; ip < num_ips; ++ip)
        {
            auto const& ip_data_vector_element = ip_data_vector[ip];
            auto const& ip_data =
                accessor_current_level_from_ip_data_vec_element(
                    ip_data_vector_element);

            if constexpr (num_comp == 1)
            {
                // scalar
                result[ip] = ip_data;
            }
            else if constexpr (num_rows == MathLib::KelvinVector::
                                               kelvin_vector_dimensions(Dim) &&
                               num_cols == 1)
            {
                // Kelvin vector
                auto const converted =
                    MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                        ip_data);

                for (unsigned comp = 0; comp < num_comp; ++comp)
                {
                    result[ip * num_comp + comp] = converted[comp];
                }
            }
            else if constexpr (num_cols == MathLib::KelvinVector::
                                               kelvin_vector_dimensions(Dim) &&
                               num_rows == 1)
            {
                static_assert(
                    num_rows != 1 /* always false in this branch */,
                    "We support Kelvin column-vectors, but not Kelvin "
                    "row-vectors. The latter are unusual and confusion with "
                    "generic vectors might be possible.");
            }
            else if constexpr (num_rows == 1 || num_cols == 1)
            {
                // row or column vector
                for (unsigned comp = 0; comp < num_comp; ++comp)
                {
                    result[ip * num_comp + comp] = ip_data[comp];
                }
            }
            else
            {
                // matrix
                // row-major traversal
                for (unsigned row = 0; row < num_rows; ++row)
                {
                    for (unsigned col = 0; col < num_cols; ++col)
                    {
                        result[ip * num_comp + row * num_cols + col] =
                            ip_data(row, col);
                    }
                }
            }
        }
        return result;
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
         &callback]<typename Class, typename Accessor>(
            ReflectionData<Class, Accessor> const& refl_data)
        {
            using MemberRef = std::invoke_result_t<Accessor, Class const&>;
            using Member = std::remove_cvref_t<MemberRef>;

            auto accessor_member_from_ip_data_vec_element =
                [accessor_next_level = refl_data.accessor,
                 accessor_current_level_from_ip_data_vec_element](
                    auto const& ip_data_vec_element) -> Member const&
            {
                return accessor_next_level(
                    accessor_current_level_from_ip_data_vec_element(
                        ip_data_vec_element));
            };

            if constexpr (is_reflectable<Member>)
            {
                forEachReflectedFlattenedIPDataAccessor<Dim>(
                    callback, detail::reflect(std::type_identity<Member>{}),
                    accessor_ip_data_vec_in_loc_asm,
                    accessor_member_from_ip_data_vec_element);
            }
            else
            {
                static_assert(is_raw_data<Member>::value,
                              "The current member is not reflectable, so we "
                              "expect it to be raw data.");

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
        [&callback]<typename Class, typename Accessor>(
            ReflectionData<Class, Accessor> const& refl_data)
        {
            static_assert(std::is_same_v<Class, LocAsmIF>,
                          "The currently processed reflection data is not for "
                          "the given LocAsmIF but for a different class.");

            using AccessorResultRef =
                std::invoke_result_t<Accessor, Class const&>;
            using AccessorResult = std::remove_cvref_t<AccessorResultRef>;

            // AccessorResult must be a std::vector<SomeType, SomeAllocator>. We
            // check that, now.
            static_assert(boost::mp11::mp_is_list<AccessorResult>::
                              value);  // std::vector<SomeType, SomeAllocator>
                                       // is a list in the Boost MP11 sense
            static_assert(
                std::is_same_v<
                    AccessorResult,
                    boost::mp11::mp_rename<AccessorResult, std::vector>>,
                "We expect a std::vector, here.");
            // Now, we know that AccessorResult is std::vector<Member>. To be
            // more specific, AccessorResult is a std::vector<IPData> and Member
            // is IPData.
            using Member = typename AccessorResult::value_type;

            auto accessor_ip_data_vec_in_loc_asm =
                [ip_data_vector_accessor =
                     refl_data.accessor](LocAsmIF const& loc_asm) -> auto const&
            { return ip_data_vector_accessor(loc_asm); };

            if constexpr (detail::is_reflectable<Member>)
            {
                detail::forEachReflectedFlattenedIPDataAccessor<Dim>(
                    callback,
                    detail::reflect(std::type_identity<Member>{}),
                    accessor_ip_data_vec_in_loc_asm);
            }
            else
            {
                static_assert(detail::is_raw_data<Member>::value,
                              "The current member is not reflectable, so we "
                              "expect it to be raw data.");

                constexpr unsigned num_comp =
                    detail::NumberOfComponents<Member>::value;

                callback(refl_data.name, num_comp,
                         detail::getFlattenedIPDataFromLocAsm<Dim>(
                             accessor_ip_data_vec_in_loc_asm));
            }
        });
}
}  // namespace ProcessLib::Reflection
