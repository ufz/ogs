/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 4, 2020, 10:27 AM
 */

#pragma once

#include <Eigen/Eigen>
#include <vector>

#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"

namespace ProcessLib
{
template <int DisplacementDim, typename IntegrationPointData,
          typename MemberType>
std::vector<double> const& getIntegrationPointKelvinVectorData(
    std::vector<IntegrationPointData,
                Eigen::aligned_allocator<IntegrationPointData>> const& ip_data,
    MemberType member, std::vector<double>& cache)
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    auto const num_intpts = ip_data.size();

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, kelvin_vector_size, num_intpts);

    for (unsigned ip = 0; ip < num_intpts; ++ip)
    {
        auto const& kelvin_vector = ip_data[ip].*member;
        cache_mat.col(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(kelvin_vector);
    }

    return cache;
}

template <int DisplacementDim, typename IntegrationPointData,
          typename MemberType>
std::vector<double> getIntegrationPointKelvinVector(
    std::vector<IntegrationPointData,
                Eigen::aligned_allocator<IntegrationPointData>> const& ip_data,
    MemberType member)
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    auto const num_intpts = ip_data.size();

    std::vector<double> ip_kelvin_vector_values;
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, Eigen::Dynamic, kelvin_vector_size, Eigen::RowMajor>>(
        ip_kelvin_vector_values, num_intpts, kelvin_vector_size);

    for (unsigned ip = 0; ip < num_intpts; ++ip)
    {
        auto const& ip_member = ip_data[ip].*member;
        cache_mat.row(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(ip_member);
    }

    return ip_kelvin_vector_values;
}

template <int DisplacementDim, typename IntegrationPointData,
          typename MemberType>
std::size_t setIntegrationPointKelvinVector(
    double const* values,
    std::vector<IntegrationPointData,
                Eigen::aligned_allocator<IntegrationPointData>>& ip_data,
    MemberType member)
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    auto const num_intpts = ip_data.size();

    auto kelvin_vector_values =
        Eigen::Map<Eigen::Matrix<double, kelvin_vector_size, Eigen::Dynamic,
                                 Eigen::ColMajor> const>(
            values, kelvin_vector_size, num_intpts);

    for (unsigned ip = 0; ip < num_intpts; ++ip)
    {
        ip_data[ip].*member =
            MathLib::KelvinVector::symmetricTensorToKelvinVector(
                kelvin_vector_values.col(ip));
    }

    return num_intpts;
}

template <typename IntegrationPointData, typename MemberType>
std::vector<double> const& getIntegrationPointScalarData(
    std::vector<IntegrationPointData,
                Eigen::aligned_allocator<IntegrationPointData>> const& ip_data,
    MemberType member, std::vector<double>& cache)
{
    auto const num_intpts = ip_data.size();

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(cache, 1,
                                                                   num_intpts);

    for (unsigned ip = 0; ip < num_intpts; ++ip)
    {
        cache_mat[ip] = ip_data[ip].*member;
    }

    return cache;
}

template <typename IntegrationPointData, typename MemberType>
std::size_t setIntegrationPointScalarData(
    double const* values,
    std::vector<IntegrationPointData,
                Eigen::aligned_allocator<IntegrationPointData>>& ip_data,
    MemberType member)
{
    auto const num_intpts = ip_data.size();

    for (unsigned ip = 0; ip < num_intpts; ++ip)
    {
        ip_data[ip].*member = values[ip];
    }
    return num_intpts;
}

}  // namespace ProcessLib
