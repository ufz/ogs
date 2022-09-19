/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"

template <int DisplacementDim>
struct IPData
{
    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

    KV kelvin;
    double scalar;
};

template <class Dim>
struct ProcessLib_IPDataAccess : ::testing::Test
{
    static constexpr int dim = Dim::value;
    static constexpr int kv_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(dim);

    static std::vector<IPData<Dim::value>> getIPData()
    {
        using KV = typename IPData<dim>::KV;

        constexpr int off_diag_size = dim == 2 ? 1 : 3;

        constexpr std::size_t num_int_pts = 10;

        std::vector<IPData<dim>> ip_data(num_int_pts);

        for (std::size_t i = 0; i < num_int_pts; ++i)
        {
            ip_data[i].kelvin =
                KV::Constant(10. * i) + KV::LinSpaced(0., kv_size - 1.);

            // compensate Kelvin vector <-> symmetric tensor conversion
            ip_data[i].kelvin.template tail<off_diag_size>() *= std::sqrt(2.0);

            ip_data[i].scalar = 10. * num_int_pts + i;
        }

        return ip_data;
    }

    static std::vector<IPData<Dim::value>> getIPDataNaNs()
    {
        using KV = typename IPData<dim>::KV;

        constexpr std::size_t num_int_pts = 10;
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();

        std::vector<IPData<dim>> ip_data(num_int_pts);

        for (std::size_t i = 0; i < num_int_pts; ++i)
        {
            ip_data[i].kelvin = KV::Constant(nan);
            ip_data[i].scalar = nan;
        }

        return ip_data;
    }

    static std::vector<double> getScalarData()
    {
        return {100, 101, 102, 103, 104, 105, 106, 107, 108, 109};
    }

    static std::vector<double> getKVDataDefaultOrder()
    {
        if constexpr (dim == 2)
        {
            return {0, 10, 20, 30, 40, 50, 60, 70, 80, 90,   // 1st comp
                    1, 11, 21, 31, 41, 51, 61, 71, 81, 91,   // 2nd comp
                    2, 12, 22, 32, 42, 52, 62, 72, 82, 92,   // 3rd comp
                    3, 13, 23, 33, 43, 53, 63, 73, 83, 93};  // 4th comp
        }
        else if constexpr (dim == 3)
        {
            return {0, 10, 20, 30, 40, 50, 60, 70, 80, 90,  //
                    1, 11, 21, 31, 41, 51, 61, 71, 81, 91,  //
                    2, 12, 22, 32, 42, 52, 62, 72, 82, 92,  //
                    3, 13, 23, 33, 43, 53, 63, 73, 83, 93,  //
                    4, 14, 24, 34, 44, 54, 64, 74, 84, 94,  //
                    5, 15, 25, 35, 45, 55, 65, 75, 85, 95};
        };
    }

    static std::vector<double> getKVDataTransposedOrder()
    {
        if constexpr (dim == 2)
        {
            return {0,  1,  2,  3,   // 1st IP
                    10, 11, 12, 13,  // 2nd IP
                    20, 21, 22, 23,  // 3rd IP
                    30, 31, 32, 33,  // ...
                    40, 41, 42, 43,  //
                    50, 51, 52, 53,  //
                    60, 61, 62, 63,  //
                    70, 71, 72, 73,  //
                    80, 81, 82, 83,  //
                    90, 91, 92, 93};
        }
        else if constexpr (dim == 3)
        {
            return {
                0,  1,  2,  3,  4,  5,   //
                10, 11, 12, 13, 14, 15,  //
                20, 21, 22, 23, 24, 25,  //
                30, 31, 32, 33, 34, 35,  //
                40, 41, 42, 43, 44, 45,  //
                50, 51, 52, 53, 54, 55,  //
                60, 61, 62, 63, 64, 65,  //
                70, 71, 72, 73, 74, 75,  //
                80, 81, 82, 83, 84, 85,  //
                90, 91, 92, 93, 94, 95,
            };
        };
    }
};

using TestCases = ::testing::Types<std::integral_constant<int, 2>,
                                   std::integral_constant<int, 3>>;

TYPED_TEST_SUITE(ProcessLib_IPDataAccess, TestCases);

TYPED_TEST(ProcessLib_IPDataAccess, GetScalarData)
{
    constexpr int dim = TypeParam::value;

    auto const ip_data = this->getIPData();

    std::vector<double> cache;

    ProcessLib::getIntegrationPointScalarData(
        ip_data, &IPData<dim>::scalar, cache);

    ASSERT_THAT(cache,
                testing::Pointwise(testing::DoubleEq(), this->getScalarData()));
}

TYPED_TEST(ProcessLib_IPDataAccess, GetKelvinVectorDataDefaultOrder)
{
    constexpr int dim = TypeParam::value;

    auto const ip_data = this->getIPData();

    std::vector<double> cache;

    ProcessLib::getIntegrationPointKelvinVectorData<dim>(
        ip_data, &IPData<dim>::kelvin, cache);

    ASSERT_THAT(
        cache,
        testing::Pointwise(testing::DoubleEq(), this->getKVDataDefaultOrder()));
}

TYPED_TEST(ProcessLib_IPDataAccess, GetKelvinVectorDataTransposedOrder)
{
    constexpr int dim = TypeParam::value;

    auto const ip_data = this->getIPData();

    const std::vector<double> cache =
        ProcessLib::getIntegrationPointKelvinVectorData<dim>(
            ip_data,
            &IPData<dim>::kelvin);  // pretty subtle: if no cache argument is
                                    // passed, data is returned transposed

    ASSERT_THAT(cache,
                testing::Pointwise(testing::DoubleEq(),
                                   this->getKVDataTransposedOrder()));
}

TYPED_TEST(ProcessLib_IPDataAccess, SetScalarData)
{
    constexpr int dim = TypeParam::value;

    auto ip_data = this->getIPDataNaNs();

    auto const cache = this->getScalarData();

    auto const num_read = ProcessLib::setIntegrationPointScalarData(
        &cache[0], ip_data, &IPData<dim>::scalar);

    ASSERT_EQ(ip_data.size(), num_read);

    auto const ip_data_expected = this->getIPData();

    for (std::size_t i = 0; i < ip_data_expected.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(ip_data_expected[i].scalar, ip_data[i].scalar)
            << "Values at integration point " << i << " differ.";
    }
}

TYPED_TEST(ProcessLib_IPDataAccess, SetKelvinVectorData)
{
    constexpr int dim = TypeParam::value;

    auto ip_data = this->getIPDataNaNs();

    auto const cache = this->getKVDataTransposedOrder();

    auto const num_read = ProcessLib::setIntegrationPointKelvinVectorData<dim>(
        &cache[0], ip_data, &IPData<dim>::kelvin);

    ASSERT_EQ(ip_data.size(), num_read);

    auto const ip_data_expected = this->getIPData();

    for (std::size_t i = 0; i < ip_data_expected.size(); ++i)
    {
        EXPECT_THAT(
            ip_data[i].kelvin,
            testing::Pointwise(testing::DoubleEq(), ip_data_expected[i].kelvin))
            << "Values at integration point " << i << " differ.";
    }
}
