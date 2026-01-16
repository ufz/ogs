// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <complex>
#include <tuple>
#include <vector>

#include "BaseLib/DemangleTypeInfo.h"

namespace
{
struct CustomStruct
{
};

class CustomClass
{
};

template <typename T>
class TemplateClass
{
};
}  // namespace

TEST(BaseLibDemangleTypeInfo, Double)
{
    auto result = BaseLib::typeToString<double>();
    EXPECT_EQ(result, "double") << result;
}

TEST(BaseLibDemangleTypeInfo, ConstDouble)
{
    auto result = BaseLib::typeToString<const double>();
    EXPECT_EQ(result, "double") << result;
}

TEST(BaseLibDemangleTypeInfo, DoubleReference)
{
    auto result = BaseLib::typeToString<double&>();
    EXPECT_EQ(result, "double") << result;
}

TEST(BaseLibDemangleTypeInfo, DoubleConstReference)
{
    auto result = BaseLib::typeToString<double const&>();
    EXPECT_EQ(result, "double") << result;
}

TEST(BaseLibDemangleTypeInfo, DoubleArray)
{
    auto result = BaseLib::typeToString<double[10]>();
    EXPECT_TRUE(result.contains("double")) << result;
    EXPECT_TRUE(result.contains("[10]")) << result;
}

TEST(BaseLibDemangleTypeInfo, Int)
{
    auto result = BaseLib::typeToString<int>();
    EXPECT_EQ(result, "int") << result;
}

TEST(BaseLibDemangleTypeInfo, ConstInt)
{
    auto result = BaseLib::typeToString<const int>();
    EXPECT_EQ(result, "int") << result;
}

TEST(BaseLibDemangleTypeInfo, IntReference)
{
    auto result = BaseLib::typeToString<int&>();
    EXPECT_EQ(result, "int") << result;
}

TEST(BaseLibDemangleTypeInfo, IntReferenceEqualsInt)
{
    auto int_result = BaseLib::typeToString<int>();
    auto int_ref_result = BaseLib::typeToString<int&>();
    // Verify that reference demangling gives same base type as non-reference
    EXPECT_EQ(int_result, int_ref_result)
        << "int should demangle to same as int&";
}

TEST(BaseLibDemangleTypeInfo, ConstIntReference)
{
    auto result = BaseLib::typeToString<const int&>();
    EXPECT_EQ(result, "int") << result;
}

TEST(BaseLibDemangleTypeInfo, PointerToInt)
{
    auto result = BaseLib::typeToString<int*>();

    EXPECT_TRUE(result.contains("*") || result.contains("ptr")) << result;
}

TEST(BaseLibDemangleTypeInfo, IntArray)
{
    auto result = BaseLib::typeToString<int[5]>();
    EXPECT_TRUE(result.contains("int")) << result;
    EXPECT_TRUE(result.contains("[5]")) << result;
}

TEST(BaseLibDemangleTypeInfo, Float)
{
    auto result = BaseLib::typeToString<float>();
    EXPECT_EQ(result, "float") << result;
}

TEST(BaseLibDemangleTypeInfo, Char)
{
    auto result = BaseLib::typeToString<char>();
    EXPECT_EQ(result, "char") << result;
}

TEST(BaseLibDemangleTypeInfo, Bool)
{
    auto result = BaseLib::typeToString<bool>();
    EXPECT_EQ(result, "bool") << result;
}

TEST(BaseLibDemangleTypeInfo, StdString)
{
    auto result = BaseLib::typeToString<std::string>();

    EXPECT_TRUE(result.contains("std::")) << result;
    EXPECT_TRUE(result.contains("basic_string<")) << result;
}

TEST(BaseLibDemangleTypeInfo, StdVector)
{
    auto result = BaseLib::typeToString<std::vector<int>>();

    EXPECT_TRUE(result.contains("std::")) << result;
    EXPECT_TRUE(result.contains("vector<int")) << result;
}

TEST(BaseLibDemangleTypeInfo, StdComplex)
{
    auto result = BaseLib::typeToString<std::complex<double>>();

    EXPECT_TRUE(result.contains("complex")) << result;
    EXPECT_TRUE(result.contains("double")) << result;
}

TEST(BaseLibDemangleTypeInfo, CustomStruct)
{
    auto result = BaseLib::typeToString<CustomStruct>();

    EXPECT_TRUE(result.contains("CustomStruct")) << result;
}

TEST(BaseLibDemangleTypeInfo, CustomClass)
{
    auto result = BaseLib::typeToString<CustomClass>();

    EXPECT_TRUE(result.contains("CustomClass")) << result;
}

TEST(BaseLibDemangleTypeInfo, TemplateClass)
{
    auto result = BaseLib::typeToString<TemplateClass<int>>();

    EXPECT_TRUE(result.contains("TemplateClass<int>")) << result;
}

TEST(BaseLibDemangleTypeInfo, StdTuple)
{
    auto result = BaseLib::typeToString<std::tuple<int, double, std::string>>();

    EXPECT_TRUE(result.contains("std::")) << result;
    EXPECT_TRUE(result.contains("tuple<")) << result;
    EXPECT_TRUE(result.contains("<int,")) << result;
    EXPECT_TRUE(result.contains("double,")) << result;
    EXPECT_TRUE(result.contains("basic_string<")) << result;
}

TEST(BaseLibDemangleTypeInfo, TupleWithReference)
{
    auto result = BaseLib::typeToString<std::tuple<int&, double>>();
    EXPECT_TRUE(result.contains("std::")) << result;
    EXPECT_TRUE(result.contains("tuple")) << result;
    EXPECT_TRUE(result.contains("int")) << result;
    EXPECT_TRUE(result.contains("&")) << result;
    EXPECT_TRUE(result.contains("double")) << result;
}

TEST(BaseLibDemangleTypeInfo, EigenVectorXd)
{
    auto result = BaseLib::typeToString<Eigen::VectorXd>();

    EXPECT_TRUE(result.contains("Eigen")) << result;
    EXPECT_TRUE(result.contains("Matrix")) << result;
    EXPECT_TRUE(result.contains("double")) << result;
    EXPECT_TRUE(result.contains("-1")) << result;
}

TEST(BaseLibDemangleTypeInfo, EigenMatrix3d)
{
    auto result = BaseLib::typeToString<Eigen::Matrix3d>();
    EXPECT_TRUE(result.contains("Eigen")) << result;
    EXPECT_TRUE(result.contains("Matrix")) << result;
    EXPECT_TRUE(result.contains("double")) << result;
    EXPECT_TRUE(result.contains("3")) << result;
}

TEST(BaseLibDemangleTypeInfo, FunctionPointer)
{
    auto result = BaseLib::typeToString<char (*)(double, int)>();
    EXPECT_TRUE(result.contains("char")) << result;
    EXPECT_TRUE(result.contains("*")) << result;
    EXPECT_TRUE(result.contains("double")) << result;
    EXPECT_TRUE(result.contains("int")) << result;
}
