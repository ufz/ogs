/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 4, 2020, 10:27 AM
 */

#pragma once

#include <Eigen/Core>
#include <span>
#include <vector>

#include "MathLib/KelvinVector.h"
#include "TransposeInPlace.h"

namespace ProcessLib
{

template <int DisplacementDim, typename IntegrationPointDataVector,
          typename IpData, typename MemberType>
std::vector<double> const& getIntegrationPointDimMatrixData(
    IntegrationPointDataVector const& ip_data_vector,
    MemberType IpData::*const member, std::vector<double>& cache)
{
    return getIntegrationPointDimMatrixData<DisplacementDim>(
        ip_data_vector,
        [member](IpData const& ip_data) -> auto const&
        { return ip_data.*member; },
        cache);
}

template <int DisplacementDim, typename IntegrationPointDataVector,
          typename Accessor>
std::vector<double> const& getIntegrationPointDimMatrixData(
    IntegrationPointDataVector const& ip_data_vector, Accessor&& accessor,
    std::vector<double>& cache)
{
    using AccessorResult = decltype(std::declval<Accessor>()(
        std::declval<IntegrationPointDataVector>()[0]));
    static_assert(std::is_lvalue_reference_v<AccessorResult>,
                  "The ip data accessor should return a reference. This check "
                  "prevents accidental copies.");

    auto const n_integration_points = ip_data_vector.size();
    constexpr int size = DisplacementDim * DisplacementDim;

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, size, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, size, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& dim_matix = accessor(ip_data_vector[ip]);
        cache_mat.col(ip) = dim_matix.template reshaped<Eigen::RowMajor>();
    }

    return cache;
}

template <int DisplacementDim, typename IntegrationPointDataVector,
          typename IpData, typename MemberType>
std::vector<double> const& getIntegrationPointKelvinVectorData(
    IntegrationPointDataVector const& ip_data_vector,
    MemberType IpData::*const member, std::vector<double>& cache)
{
    return getIntegrationPointKelvinVectorData<DisplacementDim>(
        ip_data_vector,
        [member](IpData const& ip_data) -> auto const&
        { return ip_data.*member; },
        cache);
}

template <int DisplacementDim, typename IntegrationPointDataVector,
          typename Accessor>
std::vector<double> const& getIntegrationPointKelvinVectorData(
    IntegrationPointDataVector const& ip_data_vector, Accessor&& accessor,
    std::vector<double>& cache)
{
    using AccessorResult = decltype(std::declval<Accessor>()(
        std::declval<IntegrationPointDataVector>()[0]));
    static_assert(std::is_lvalue_reference_v<AccessorResult>,
                  "The ip data accessor should return a reference. This check "
                  "prevents accidental copies.");

    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    auto const n_integration_points = ip_data_vector.size();

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, kelvin_vector_size, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& kelvin_vector = accessor(ip_data_vector[ip]);
        cache_mat.col(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(kelvin_vector);
    }

    return cache;
}

//! Overload without \c cache argument.
//!
//! \note This function returns the data in transposed storage order compared to
//! the overloads that have a \c cache argument.
template <int DisplacementDim, typename IntegrationPointDataVector,
          typename MemberType>
std::vector<double> getIntegrationPointKelvinVectorData(
    IntegrationPointDataVector const& ip_data_vector, MemberType member)
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

    return transposeInPlace<kelvin_vector_size>(
        [&](std::vector<double>& values)
        {
            return getIntegrationPointKelvinVectorData<DisplacementDim>(
                ip_data_vector, member, values);
            ;
        });
}

template <int DisplacementDim, typename IntegrationPointDataVector,
          typename IpData, typename MemberType>
std::size_t setIntegrationPointKelvinVectorData(
    double const* values,
    IntegrationPointDataVector& ip_data_vector,
    MemberType IpData::*const member)
{
    return setIntegrationPointKelvinVectorData<DisplacementDim>(
        values, ip_data_vector,
        [member](IpData& ip_data) -> auto& { return ip_data.*member; });
}

template <int DisplacementDim, typename IntegrationPointDataVector,
          typename Accessor>
std::size_t setIntegrationPointKelvinVectorData(
    double const* values,
    IntegrationPointDataVector& ip_data_vector,
    Accessor&& accessor)
{
    using AccessorResult = decltype(std::declval<Accessor>()(
        std::declval<IntegrationPointDataVector>()[0]));
    static_assert(std::is_lvalue_reference_v<AccessorResult>,
                  "The ip data accessor must return a reference.");

    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    auto const n_integration_points = ip_data_vector.size();

    auto kelvin_vector_values =
        Eigen::Map<Eigen::Matrix<double, kelvin_vector_size, Eigen::Dynamic,
                                 Eigen::ColMajor> const>(
            values, kelvin_vector_size, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        accessor(ip_data_vector[ip]) =
            MathLib::KelvinVector::symmetricTensorToKelvinVector(
                kelvin_vector_values.col(ip));
    }

    return n_integration_points;
}

template <typename IntegrationPointDataVector, typename IpData,
          typename MemberType>
std::vector<double> const& getIntegrationPointScalarData(
    IntegrationPointDataVector const& ip_data_vector,
    MemberType IpData::*const member, std::vector<double>& cache)
{
    return getIntegrationPointScalarData(
        ip_data_vector,
        [member](IpData const& ip_data) -> auto const&
        { return ip_data.*member; },
        cache);
}

template <typename IntegrationPointDataVector, typename Accessor>
std::vector<double> const& getIntegrationPointScalarData(
    IntegrationPointDataVector const& ip_data_vector, Accessor&& accessor,
    std::vector<double>& cache)
{
    using AccessorResult = decltype(std::declval<Accessor>()(
        std::declval<IntegrationPointDataVector>()[0]));
    static_assert(std::is_lvalue_reference_v<AccessorResult>,
                  "The ip data accessor should return a reference. This check "
                  "prevents accidental copies.");

    auto const n_integration_points = ip_data_vector.size();

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, 1, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        cache_mat[ip] = accessor(ip_data_vector[ip]);
    }

    return cache;
}

template <typename IntegrationPointDataVector, typename IpData,
          typename MemberType>
std::size_t setIntegrationPointScalarData(
    double const* values,
    IntegrationPointDataVector& ip_data_vector,
    MemberType IpData::*const member)
{
    return setIntegrationPointScalarData(values, ip_data_vector,
                                         [member](IpData& ip_data) -> auto&
                                         { return ip_data.*member; });
}

template <typename IntegrationPointDataVector, typename Accessor>
std::size_t setIntegrationPointScalarData(
    double const* values,
    IntegrationPointDataVector& ip_data_vector,
    Accessor&& accessor)
{
    using AccessorResult = decltype(std::declval<Accessor>()(
        std::declval<IntegrationPointDataVector>()[0]));
    static_assert(std::is_lvalue_reference_v<AccessorResult>,
                  "The ip data accessor must return a reference.");

    auto const n_integration_points = ip_data_vector.size();

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        accessor(ip_data_vector[ip]) = values[ip];
    }
    return n_integration_points;
}

template <typename IntegrationPointDataVector, typename MemberType,
          typename MaterialStateVariables>
std::vector<double> getIntegrationPointDataMaterialStateVariables(
    IntegrationPointDataVector const& ip_data_vector,
    MemberType member,
    std::function<std::span<double>(MaterialStateVariables&)>
        get_values_span,
    int const n_components)
{
    std::vector<double> result;
    result.reserve(ip_data_vector.size() * n_components);

    for (auto& ip_data : ip_data_vector)
    {
        auto const values_span = get_values_span(*(ip_data.*member));
        assert(values_span.size() == static_cast<std::size_t>(n_components));

        result.insert(end(result), values_span.begin(), values_span.end());
    }

    return result;
}

template <typename IntegrationPointDataVector, typename MemberType,
          typename MaterialStateVariables>
std::size_t setIntegrationPointDataMaterialStateVariables(
    double const* values,
    IntegrationPointDataVector& ip_data_vector,
    MemberType member,
    std::function<std::span<double>(MaterialStateVariables&)>
        get_values_span)
{
    auto const n_integration_points = ip_data_vector.size();

    std::size_t position = 0;
    for (auto const& ip_data : ip_data_vector)
    {
        auto const values_span = get_values_span(*(ip_data.*member));
        std::copy_n(values + position, values_span.size(), values_span.begin());
        position += values_span.size();
    }
    return n_integration_points;
}
}  // namespace ProcessLib
