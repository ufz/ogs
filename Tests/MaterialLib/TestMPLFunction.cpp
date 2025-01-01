/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <limits>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Properties/Function.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/VectorizedTensor.h"

namespace MPL = MaterialPropertyLib;

struct MPLFunction : public ::testing::Test
{
    static double constexpr nan = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;

    using KV2 = MathLib::KelvinVector::KelvinVectorType<2>;
    using KV3 = MathLib::KelvinVector::KelvinVectorType<3>;

    using VT2 = MathLib::VectorizedTensor::Type<2>;
    using VT3 = MathLib::VectorizedTensor::Type<3>;
};

TEST_F(MPLFunction, ScalarScalar)
{
    vars.temperature = 2.;

    MPL::Property const& f = MPL::Function{
        "test_function", {"temperature"}, {{"temperature", {"1"}}}};
    ASSERT_EQ(2., f.value<double>(vars, {}, nan, nan));
    ASSERT_EQ(1.,
              f.dValue<double>(vars, MPL::Variable::temperature, {}, nan, nan));
}

TEST_F(MPLFunction, ScalarVector)
{
    vars.temperature = 2.;

    MPL::Property const& f =
        MPL::Function{"test_function",
                      {"temperature", "temperature^2"},
                      {{"temperature", {"1", "2*temperature"}}}};
    ASSERT_EQ((Eigen::Vector2d{2., 4.}),
              (f.value<Eigen::Vector2d>(vars, {}, nan, nan)));
    ASSERT_EQ((Eigen::Vector2d{1., 4.}),
              f.dValue<Eigen::Vector2d>(
                  vars, MPL::Variable::temperature, {}, nan, nan));
}

TEST_F(MPLFunction, KelvinVector2Scalar)
{
    vars.stress.emplace<KV2>(1, 2, 3, 4 * std::sqrt(2.));
    vars.temperature = 273.15;

    MPL::Property const& f =
        MPL::Function{"test_function", {"avg(stress) * temperature"}, {}};
    ASSERT_EQ(((1 + 4) * 0.5) * 273.15, f.value<double>(vars, {}, nan, nan));
}

TEST_F(MPLFunction, KelvinVector3Scalar)
{
    vars.stress.emplace<KV3>(
        1, 2, 3, 4 * std::sqrt(2.), 5 * std::sqrt(2.), 6 * std::sqrt(2.));
    vars.temperature = 273.15;

    MPL::Property const& f =
        MPL::Function{"test_function", {"avg(stress) * temperature"}, {}};
    ASSERT_EQ(((1 + 6) * 0.5) * 273.15, f.value<double>(vars, {}, nan, nan));
}

TEST_F(MPLFunction, KelvinVector23Scalar)
{
    // Mixed dimensions of vectorial quantities are expected to fail.
    vars.stress.emplace<KV3>(1, 2, 3, 4, 5, 6);
    vars.total_strain.emplace<KV2>(1, 2, 3, 4);

    MPL::Property const& f =
        MPL::Function{"test_function", {"avg(stress) * avg(total_strain)"}, {}};
    ASSERT_ANY_THROW(f.value<double>(vars, {}, nan, nan));
}

TEST_F(MPLFunction, VectorizedTensor2Scalar)
{
    vars.deformation_gradient.emplace<VT2>(1, 2, 3, 4, 5);
    vars.temperature = 273.15;

    MPL::Property const& f = MPL::Function{
        "test_function", {"avg(deformation_gradient) * temperature"}, {}};
    ASSERT_EQ(((1 + 5) * 0.5) * 273.15, f.value<double>(vars, {}, nan, nan));
}

TEST_F(MPLFunction, VectorizedTensor3Scalar)
{
    vars.deformation_gradient.emplace<VT3>(1, 2, 3, 4, 5, 6, 7, 8, 9);
    vars.temperature = 273.15;

    MPL::Property const& f = MPL::Function{
        "test_function", {"avg(deformation_gradient) * temperature"}, {}};
    ASSERT_EQ(((1 + 9) * 0.5) * 273.15, f.value<double>(vars, {}, nan, nan));
}

TEST_F(MPLFunction, ScalarUninitialized)
{
    // The vars.temperature = is not initialized.

    MPL::Property const& f = MPL::Function{
        "test_function", {"temperature"}, {{"temperature", {"1"}}}};
    ASSERT_ANY_THROW(f.value<double>(vars, {}, nan, nan));
    ASSERT_ANY_THROW(
        f.dValue<double>(vars, MPL::Variable::temperature, {}, nan, nan));
}

TEST_F(MPLFunction, KelvinVector2Uninitialized)
{
    vars.temperature = 273.15;
    // The vars.stress is not initialized.

    MPL::Property const& f =
        MPL::Function{"test_function", {"avg(stress) * temperature"}, {}};
    ASSERT_ANY_THROW(f.value<double>(vars, {}, nan, nan));
}

TEST_F(MPLFunction, VectorizedTensor2Uninitialized)
{
    vars.temperature = 273.15;
    // The vars.deformation_gradient is not initialized.

    MPL::Property const& f = MPL::Function{
        "test_function", {"avg(deformation_gradient) * temperature"}, {}};
    ASSERT_ANY_THROW(f.value<double>(vars, {}, nan, nan));
}

TEST_F(MPLFunction, KelvinVector3Uninitialized)
{
    vars.temperature = 273.15;
    // The vars.stress is not initialized.

    MPL::Property const& f =
        MPL::Function{"test_function", {"avg(stress) * temperature"}, {}};
    ASSERT_ANY_THROW(f.value<double>(vars, {}, nan, nan));
}

TEST_F(MPLFunction, VectorizedTensor3Uninitialized)
{
    vars.temperature = 273.15;
    // The vars.deformation_gradient is not initialized.

    MPL::Property const& f = MPL::Function{
        "test_function", {"avg(deformation_gradient) * temperature"}, {}};
    ASSERT_ANY_THROW(f.value<double>(vars, {}, nan, nan));
}
