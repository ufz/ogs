/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <Eigen/Core>

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

using ProcessLib_IPDataAccess_TestCases =
    ::testing::Types<std::integral_constant<int, 2>,
                     std::integral_constant<int, 3>>;

TYPED_TEST_SUITE(ProcessLib_IPDataAccess, ProcessLib_IPDataAccess_TestCases);

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

template <int DisplacementDim>
struct IPDimMatrixData
{
    Eigen::Matrix<double, DisplacementDim, DisplacementDim, Eigen::RowMajor>
        dim_matrix_row_major;
    Eigen::Matrix<double, DisplacementDim, DisplacementDim, Eigen::ColMajor>
        dim_matrix_col_major;
};

template <class Dim>
struct ProcessLib_IPDimMatrixDataAccess : ::testing::Test
{
    static constexpr int dim = Dim::value;

    static std::vector<IPDimMatrixData<Dim::value>> getIPData()
    {
        constexpr std::size_t num_int_pts = 4;

        std::vector<IPDimMatrixData<dim>> ip_data(num_int_pts);

        for (std::size_t i = 0; i < num_int_pts; ++i)
        {
            // Create a test square matrix from vector.
            constexpr auto N = dim * dim;
            constexpr auto stride = 10;  // The data for each IP will be 10
                                         // apart from the previous/next IP.
            static_assert(stride >= N);

            double const low = static_cast<double>(stride * i);
            double const high = static_cast<double>(low + N - 1);
            auto const K0 = Eigen::VectorXd::LinSpaced(N, low, high)
                                .reshaped<Eigen::RowMajor>(dim, dim)
                                .eval();

            // Double-check row major storage order.
            if constexpr (dim > 1)
            {
                EXPECT_DOUBLE_EQ(low + 1, K0(0, 1))
                    << "Something extremely bad happened: there is an error in "
                       "the internal logic of this test case.";
            }

            ip_data[i].dim_matrix_row_major = K0;
            ip_data[i].dim_matrix_col_major = K0;
        }

        return ip_data;
    }

    static std::vector<double> getDimMatrixData()
    {
        if constexpr (dim == 1)
        {
            return {0, 10, 20, 30};  // K(0,0), ip = 0, 1, 2, 3
        }
        if constexpr (dim == 2)
        {
            return {
                0, 10, 20, 30,  // K(0,0), ip = 0, 1, 2, 3
                1, 11, 21, 31,  // K(0,1), ip = 0, 1, 2, 3
                2, 12, 22, 32,  // K(1,0), ip = 0, 1, 2, 3
                3, 13, 23, 33   // K(1,1), ip = 0, 1, 2, 3
            };
        }
        else if constexpr (dim == 3)
        {
            return {0, 10, 20, 30,  // K(0,0), ip = 0, 1, 2, 3
                    1, 11, 21, 31,  // K(0,1), ip = 0, 1, 2, 3
                    2, 12, 22, 32,  // K(0,2), ip = 0, 1, 2, 3
                    3, 13, 23, 33,  // K(1,0), ip = 0, 1, 2, 3
                    4, 14, 24, 34,  // ...
                    5, 15, 25, 35,  //
                    6, 16, 26, 36,  //
                    7, 17, 27, 37,  //
                    8, 18, 28, 38};
        }
    }
};

using ProcessLib_IPDimMatrixDataAccess_TestCases =
    ::testing::Types<std::integral_constant<int, 1>,
                     std::integral_constant<int, 2>,
                     std::integral_constant<int, 3>>;

TYPED_TEST_SUITE(ProcessLib_IPDimMatrixDataAccess,
                 ProcessLib_IPDimMatrixDataAccess_TestCases);

TYPED_TEST(ProcessLib_IPDimMatrixDataAccess, GetDimMatrixData)
{
    constexpr int dim = TypeParam::value;

    auto const ip_data = this->getIPData();

    std::vector<double> cache;

    ProcessLib::getIntegrationPointDimMatrixData<dim>(
        ip_data, &IPDimMatrixData<dim>::dim_matrix_row_major, cache);

    ASSERT_THAT(
        cache,
        testing::Pointwise(testing::DoubleEq(), this->getDimMatrixData()));

    ProcessLib::getIntegrationPointDimMatrixData<dim>(
        ip_data, &IPDimMatrixData<dim>::dim_matrix_col_major, cache);

    ASSERT_THAT(
        cache,
        testing::Pointwise(testing::DoubleEq(), this->getDimMatrixData()));
}
