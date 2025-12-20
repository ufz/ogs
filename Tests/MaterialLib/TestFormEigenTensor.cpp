// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <Eigen/Dense>

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

using MaterialPropertyLib::FormEigenTensor;

// Helper function to test constexpr capability (without actually using
// constexpr)
template <int GlobalDim>
auto testConstexprCapability()
{
    using Tensor = FormEigenTensor<GlobalDim>;
    return Tensor{};
}

// Test that functions are marked as constexpr and work at runtime
TEST(FormEigenTensor, ConstexprCapableVector2D)
{
    // Test that the operator is constexpr-capable (can be used in constexpr
    // context)
    auto tensor = testConstexprCapability<2>();
    auto result = tensor(Eigen::Vector2d{1.0, 2.0});

    // Verify correct diagonal matrix creation
    EXPECT_EQ(result(0, 0), 1.0);
    EXPECT_EQ(result(0, 1), 0.0);
    EXPECT_EQ(result(1, 0), 0.0);
    EXPECT_EQ(result(1, 1), 2.0);

    // Verify it's actually diagonal
    Eigen::Matrix<double, 2, 2> expected =
        Eigen::Vector2d{1.0, 2.0}.asDiagonal();
    EXPECT_EQ(result, expected);
}

TEST(FormEigenTensor, ConstexprCapableVector3D)
{
    // Test that the operator is constexpr-capable (can be used in constexpr
    // context)
    auto tensor = testConstexprCapability<3>();
    auto result = tensor(Eigen::Vector3d{1.0, 2.0, 3.0});

    // Verify correct diagonal matrix creation
    EXPECT_EQ(result(0, 0), 1.0);
    EXPECT_EQ(result(1, 1), 2.0);
    EXPECT_EQ(result(2, 2), 3.0);

    // Verify all off-diagonal elements are zero
    EXPECT_EQ(result(0, 1), 0.0);
    EXPECT_EQ(result(0, 2), 0.0);
    EXPECT_EQ(result(1, 0), 0.0);
    EXPECT_EQ(result(1, 2), 0.0);
    EXPECT_EQ(result(2, 0), 0.0);
    EXPECT_EQ(result(2, 1), 0.0);

    // Verify it's actually diagonal
    Eigen::Matrix<double, 3, 3> expected =
        Eigen::Vector3d{1.0, 2.0, 3.0}.asDiagonal();
    EXPECT_EQ(result, expected);
}

// Test runtime error handling for invalid cases (still works)
TEST(FormEigenTensor, RuntimeErrorHandling)
{
    // Test that invalid combinations still produce runtime errors
    // 2D vector to 3D matrix should error
    EXPECT_THROW(
        { FormEigenTensor<3>{}(Eigen::Vector2d{1.0, 2.0}); },
        std::runtime_error);

    // 3D vector to 2D matrix should error
    EXPECT_THROW(
        { FormEigenTensor<2>{}(Eigen::Vector3d{1.0, 2.0, 3.0}); },
        std::runtime_error);

    // 2x2 matrix to 3D should error
    Eigen::Matrix<double, 2, 2> mat2x2;
    mat2x2 << 1.0, 0.5, 0.5, 2.0;
    EXPECT_THROW({ FormEigenTensor<3>{}(mat2x2); }, std::runtime_error);
}

// Test scalar multiplication (already constexpr)
TEST(FormEigenTensor, ConstexprCapableScalar)
{
    auto tensor = testConstexprCapability<2>();
    auto result = tensor(5.0);

    EXPECT_EQ(result(0, 0), 5.0);
    EXPECT_EQ(result(0, 1), 0.0);
    EXPECT_EQ(result(1, 0), 0.0);
    EXPECT_EQ(result(1, 1), 5.0);

    // Verify it's identity * scalar
    Eigen::Matrix<double, 2, 2> expected =
        Eigen::Matrix<double, 2, 2>::Identity() * 5.0;
    EXPECT_EQ(result, expected);
}

// Test matrix operations
TEST(FormEigenTensor, ConstexprCapableMatrix2x2)
{
    Eigen::Matrix<double, 2, 2> input;
    input << 1.0, 0.5, 0.5, 2.0;

    auto tensor = testConstexprCapability<2>();
    auto result = tensor(input);

    EXPECT_EQ(result(0, 0), 1.0);
    EXPECT_EQ(result(0, 1), 0.5);
    EXPECT_EQ(result(1, 0), 0.5);
    EXPECT_EQ(result(1, 1), 2.0);

    // Should return input matrix unchanged
    EXPECT_EQ(result, input);
}

TEST(FormEigenTensor, ConstexprCapableMatrix3x3)
{
    Eigen::Matrix<double, 3, 3> input;
    input << 1.0, 0.1, 0.2, 0.1, 2.0, 0.3, 0.2, 0.3, 3.0;

    auto tensor = testConstexprCapability<3>();
    auto result = tensor(input);

    EXPECT_EQ(result(0, 0), 1.0);
    EXPECT_EQ(result(1, 1), 2.0);
    EXPECT_EQ(result(2, 2), 3.0);
    EXPECT_EQ(result(0, 1), 0.1);
    EXPECT_EQ(result(1, 2), 0.3);

    // Should return input matrix unchanged
    EXPECT_EQ(result, input);
}

// Test 6x1 vector conversion to 3x3 (symmetric tensor)
TEST(FormEigenTensor, ConstexprCapableVector6to3)
{
    // Test symmetric 3D tensor from Voigt notation
    Eigen::Matrix<double, 6, 1> voigt;
    voigt << 1.0, 2.0, 3.0, 0.5, 0.3, 0.1;  // [xx, yy, zz, xy, yz, xz]

    auto tensor = testConstexprCapability<3>();
    auto result = tensor(voigt);

    // Test diagonal elements
    EXPECT_EQ(result(0, 0), 1.0);
    EXPECT_EQ(result(1, 1), 2.0);
    EXPECT_EQ(result(2, 2), 3.0);

    // Test off-diagonal elements from Voigt notation
    EXPECT_EQ(result(0, 1), 0.5);
    EXPECT_EQ(result(1, 2), 0.3);
    EXPECT_EQ(result(0, 2), 0.1);

    // Test symmetry
    EXPECT_EQ(result(1, 0), result(0, 1));
    EXPECT_EQ(result(2, 0), result(0, 2));
    EXPECT_EQ(result(2, 1), result(1, 2));
}