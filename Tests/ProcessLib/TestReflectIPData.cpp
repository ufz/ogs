/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <numeric>
#include <random>

#include "ProcessLib/Reflection/ReflectionIPData.h"
#include "ProcessLib/Reflection/ReflectionSetIPData.h"

template <int Dim>
struct Level3
{
    MathLib::KelvinVector::KelvinVectorType<Dim> kelvin3;
    Eigen::Vector<double, Dim> vector3;
    double scalar3;
    Eigen::Matrix<double, Dim, 4, Eigen::RowMajor> matrix3;
    Eigen::Matrix<double, 4, Dim, Eigen::RowMajor> matrix3_1;
    // Same number of components as Kelvin vector in 2D. Test that Kelvin vector
    // and matrix code are not mixed up.
    Eigen::Matrix<double, 2, 2, Eigen::RowMajor> matrix3_2;

    static auto reflect()
    {
        using namespace ProcessLib::Reflection;
        return std::tuple{makeReflectionData("kelvin3", &Level3::kelvin3),
                          makeReflectionData("vector3", &Level3::vector3),
                          makeReflectionData("scalar3", &Level3::scalar3),
                          makeReflectionData("matrix3", &Level3::matrix3),
                          makeReflectionData("matrix3_1", &Level3::matrix3_1),
                          makeReflectionData("matrix3_2", &Level3::matrix3_2)};
    }
};

template <int Dim>
struct Level3b
{
    double scalar3b;

    static auto reflect()
    {
        using namespace ProcessLib::Reflection;
        return std::tuple{makeReflectionData("scalar3b", &Level3b::scalar3b)};
    }
};

template <int Dim>
using Level2 = std::tuple<Level3<Dim>, Level3b<Dim>>;

template <int Dim>
struct Level2b
{
    double scalar2b;

    static auto reflect()
    {
        using namespace ProcessLib::Reflection;
        return std::tuple{makeReflectionData("scalar2b", &Level2b::scalar2b)};
    }
};

template <int Dim>
struct Level1
{
    MathLib::KelvinVector::KelvinVectorType<Dim> kelvin1;
    Eigen::Vector<double, Dim> vector1;
    double scalar1;
    Level2<Dim> level2;
    Level2b<Dim> level2b;

    static auto reflect()
    {
        using namespace ProcessLib::Reflection;
        return std::tuple{makeReflectionData("kelvin1", &Level1::kelvin1),
                          makeReflectionData("vector1", &Level1::vector1),
                          makeReflectionData("scalar1", &Level1::scalar1),
                          makeReflectionData(&Level1::level2),
                          makeReflectionData(&Level1::level2b)};
    }
};

template <int Dim>
struct Level1b
{
    double scalar1b;

    static auto reflect()
    {
        using namespace ProcessLib::Reflection;
        return std::tuple{makeReflectionData("scalar1b", &Level1b::scalar1b)};
    }
};

template <int Dim>
struct LocAsmIF
{
    explicit LocAsmIF(unsigned const num_ips)
        : ip_data_scalar(num_ips),
          ip_data_vector(num_ips),
          ip_data_kelvin(num_ips),
          ip_data_level1(num_ips),
          ip_data_level1b(num_ips)
    {
    }

    std::size_t numIPs() const { return ip_data_scalar.size(); }

    std::vector<double> ip_data_scalar;
    std::vector<Eigen::Vector<double, Dim>> ip_data_vector;
    std::vector<MathLib::KelvinVector::KelvinVectorType<Dim>> ip_data_kelvin;
    std::vector<Level1<Dim>> ip_data_level1;
    std::vector<Level1b<Dim>> ip_data_level1b;

    static auto reflect()
    {
        using namespace ProcessLib::Reflection;
        return std::tuple{
            makeReflectionData("scalar", &LocAsmIF::ip_data_scalar),
            makeReflectionData("vector", &LocAsmIF::ip_data_vector),
            makeReflectionData("kelvin", &LocAsmIF::ip_data_kelvin),
            makeReflectionData(&LocAsmIF::ip_data_level1),
            makeReflectionData(&LocAsmIF::ip_data_level1b)};
    }
};

template <int dim>
struct NumCompAndFunction
{
    unsigned num_comp;
    std::function<std::vector<double>(LocAsmIF<dim> const&)> function;
};

// Prepares scalar IP data for the passed local assembler.
//
// The IP data are a sequence of double values starting at the passed start
// value and incremented by one for each integration point.
//
// The location of the prepared data is specified by the IP data accessor
// callback function.
//
// Returns the expected data for use in unit test checks.
template <int dim>
std::vector<double> initScalar(LocAsmIF<dim>& loc_asm,
                               double const start_value,
                               auto const ip_data_accessor,
                               bool const for_read_test)
{
    auto const num_int_pts = loc_asm.numIPs();

    // init ip data in the local assembler
    if (for_read_test)
    {
        for (std::size_t ip = 0; ip < num_int_pts; ++ip)
        {
            ip_data_accessor(loc_asm, ip) = start_value + ip;
        }
    }
    else
    {
        for (std::size_t ip = 0; ip < num_int_pts; ++ip)
        {
            ip_data_accessor(loc_asm, ip) =
                std::numeric_limits<double>::quiet_NaN();
        }
    }

    // prepare reference data
    std::vector<double> scalar_expected(num_int_pts);
    iota(begin(scalar_expected), end(scalar_expected), start_value);
    return scalar_expected;
}

// Prepares vector valued IP data for the passed local assembler.
//
// The IP data are a sequence of double values starting at the passed start
// value and incremented by one for each integration point and vector
// component.
//
// The location of the prepared data is specified by the IP data  accessor
// callback function.
//
// Returns the expected data for use in unit test checks.
template <int dim>
std::vector<double> initVector(LocAsmIF<dim>& loc_asm,
                               double const start_value,
                               auto const ip_data_accessor,
                               bool const for_read_test)
{
    auto const num_int_pts = loc_asm.numIPs();

    // init ip data in the local assembler
    if (for_read_test)
    {
        for (std::size_t ip = 0; ip < num_int_pts; ++ip)
        {
            ip_data_accessor(loc_asm, ip) =
                Eigen::Vector<double, dim>::LinSpaced(
                    dim, ip * dim + start_value,
                    ip * dim + start_value - 1 + dim);
        }
    }
    else
    {
        for (std::size_t ip = 0; ip < num_int_pts; ++ip)
        {
            ip_data_accessor(loc_asm, ip) =
                Eigen::Vector<double, dim>::Constant(
                    std::numeric_limits<double>::quiet_NaN());
        }
    }

    // prepare reference data
    std::vector<double> vector_expected(num_int_pts * dim);
    iota(begin(vector_expected), end(vector_expected), start_value);
    return vector_expected;
}

// Prepares Kelvin vector valued IP data for the passed local assembler.
//
// The IP data are a sequence of double values starting at the passed start
// value and incremented by one for each integration point and Kelvin vector
// component.
//
// The location of the prepared data is specified by the IP data  accessor
// callback function.
//
// Returns the expected data for use in unit test checks.
template <int dim>
std::vector<double> initKelvin(LocAsmIF<dim>& loc_asm,
                               double const start_value,
                               auto const ip_data_accessor,
                               bool const for_read_test)
{
    auto constexpr kv_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(dim);

    auto const num_int_pts = loc_asm.numIPs();

    // init ip data in the local assembler
    if (for_read_test)
    {
        for (std::size_t ip = 0; ip < num_int_pts; ++ip)
        {
            ip_data_accessor(loc_asm, ip) =
                MathLib::KelvinVector::symmetricTensorToKelvinVector(
                    Eigen::Vector<double, kv_size>::LinSpaced(
                        kv_size, ip * kv_size + start_value,
                        ip * kv_size + start_value - 1 + kv_size));
        }
    }
    else
    {
        for (std::size_t ip = 0; ip < num_int_pts; ++ip)
        {
            ip_data_accessor(loc_asm, ip) =
                Eigen::Vector<double, kv_size>::Constant(
                    std::numeric_limits<double>::quiet_NaN());
        }
    }

    // prepare reference data
    std::vector<double> vector_expected(num_int_pts * kv_size);
    iota(begin(vector_expected), end(vector_expected), start_value);
    return vector_expected;
}

// Prepares matrix valued IP data for the passed local assembler.
//
// The IP data are a sequence of double values starting at the passed start
// value and incremented by one for each integration point and matrix entry.
//
// The location of the prepared data is specified by the IP data  accessor
// callback function.
//
// Returns the expected data for use in unit test checks.
template <int dim>
std::vector<double> initMatrix(LocAsmIF<dim>& loc_asm,
                               double const start_value,
                               auto const ip_data_accessor,
                               bool const for_read_test)
{
    using MatrixType = std::remove_cvref_t<
        std::invoke_result_t<std::remove_cvref_t<decltype(ip_data_accessor)>,
                             LocAsmIF<dim> const&, unsigned /* ip */>>;
    auto constexpr rows = MatrixType::RowsAtCompileTime;
    auto constexpr cols = MatrixType::ColsAtCompileTime;
    auto constexpr size = rows * cols;

    auto const num_int_pts = loc_asm.numIPs();

    // init ip data in the local assembler
    if (for_read_test)
    {
        for (std::size_t ip = 0; ip < num_int_pts; ++ip)
        {
            ip_data_accessor(loc_asm, ip) =
                Eigen::Vector<double, size>::LinSpaced(
                    size, ip * size + start_value,
                    ip * size + start_value - 1 + size)
                    .template reshaped<Eigen::RowMajor>(rows, cols);
        }
    }
    else
    {
        for (std::size_t ip = 0; ip < num_int_pts; ++ip)
        {
            ip_data_accessor(loc_asm, ip) =
                Eigen::Vector<double, size>::Constant(
                    std::numeric_limits<double>::quiet_NaN())
                    .template reshaped<Eigen::RowMajor>(rows, cols);
        }
    }

    // prepare reference data
    std::vector<double> vector_expected(num_int_pts * size);
    iota(begin(vector_expected), end(vector_expected), start_value);
    return vector_expected;
}

template <int dim>
struct ReferenceData
{
    std::vector<double> scalar;
    std::vector<double> vector;
    std::vector<double> kelvin;

    std::vector<double> scalar1;
    std::vector<double> vector1;
    std::vector<double> kelvin1;

    std::vector<double> scalar3;
    std::vector<double> vector3;
    std::vector<double> kelvin3;
    std::vector<double> matrix3;
    std::vector<double> matrix3_1;
    std::vector<double> matrix3_2;

    std::vector<double> scalar1b;
    std::vector<double> scalar2b;
    std::vector<double> scalar3b;

    static ReferenceData<dim> create(LocAsmIF<dim>& loc_asm,
                                     bool const for_read_test)
    {
        std::random_device ran_dev;
        std::mt19937 ran_gen(ran_dev());
        std::uniform_real_distribution<> ran_dist(1.0, 2.0);
        auto start_value = [&]() { return ran_dist(ran_gen); };

        ReferenceData<dim> ref;

        // level 0 - data preparation //////////////////////////////////////////

        ref.scalar = initScalar(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return loc_asm.ip_data_scalar[ip];
            },
            for_read_test);

        ref.vector = initVector(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return loc_asm.ip_data_vector[ip];
            },
            for_read_test);

        ref.kelvin = initKelvin(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return loc_asm.ip_data_kelvin[ip];
            },
            for_read_test);

        // level 1 - data preparation //////////////////////////////////////////

        ref.scalar1 = initScalar(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return loc_asm.ip_data_level1[ip].scalar1;
            },
            for_read_test);

        ref.vector1 = initVector(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return loc_asm.ip_data_level1[ip].vector1;
            },
            for_read_test);

        ref.kelvin1 = initKelvin(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return loc_asm.ip_data_level1[ip].kelvin1;
            },
            for_read_test);

        // level 3 - data preparation //////////////////////////////////////////

        ref.scalar3 = initScalar(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return std::get<Level3<dim>>(loc_asm.ip_data_level1[ip].level2)
                    .scalar3;
            },
            for_read_test);

        ref.vector3 = initVector(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return std::get<Level3<dim>>(loc_asm.ip_data_level1[ip].level2)
                    .vector3;
            },
            for_read_test);

        ref.kelvin3 = initKelvin(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return std::get<Level3<dim>>(loc_asm.ip_data_level1[ip].level2)
                    .kelvin3;
            },
            for_read_test);

        ref.matrix3 = initMatrix(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return std::get<Level3<dim>>(loc_asm.ip_data_level1[ip].level2)
                    .matrix3;
            },
            for_read_test);

        ref.matrix3_1 = initMatrix(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return std::get<Level3<dim>>(loc_asm.ip_data_level1[ip].level2)
                    .matrix3_1;
            },
            for_read_test);

        ref.matrix3_2 = initMatrix(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return std::get<Level3<dim>>(loc_asm.ip_data_level1[ip].level2)
                    .matrix3_2;
            },
            for_read_test);

        // b levels - data preparation /////////////////////////////////////////
        // b levels test that the reflection implementation recurses on multiple
        // members, not only on one.

        ref.scalar1b = initScalar(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return loc_asm.ip_data_level1b[ip].scalar1b;
            },
            for_read_test);

        ref.scalar2b = initScalar(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return loc_asm.ip_data_level1[ip].level2b.scalar2b;
            },
            for_read_test);

        ref.scalar3b = initScalar(
            loc_asm, start_value(),
            [](auto& loc_asm, unsigned const ip) -> auto& {
                return std::get<Level3b<dim>>(loc_asm.ip_data_level1[ip].level2)
                    .scalar3b;
            },
            for_read_test);

        return ref;
    }
};

template <class Dim>
struct ProcessLib_ReflectIPData : ::testing::Test
{
    static constexpr auto dim = Dim::value;
};

using ProcessLib_ReflectIPData_TestCases =
    ::testing::Types<std::integral_constant<int, 2>,
                     std::integral_constant<int, 3>>;

TYPED_TEST_SUITE(ProcessLib_ReflectIPData, ProcessLib_ReflectIPData_TestCases);

TYPED_TEST(ProcessLib_ReflectIPData, ReadTest)
{
    constexpr int dim = TypeParam::value;
    auto constexpr kv_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(dim);

    using LocAsm = LocAsmIF<dim>;

    std::size_t const num_int_pts = 8;
    LocAsm loc_asm(num_int_pts);

    auto const ref = ReferenceData<dim>::create(loc_asm, true);

    // function under test /////////////////////////////////////////////////////

    std::map<std::string, NumCompAndFunction<dim>>
        map_name_to_num_comp_and_function;

    ProcessLib::Reflection::forEachReflectedFlattenedIPDataAccessor<dim,
                                                                    LocAsm>(
        LocAsm::reflect(),
        [&map_name_to_num_comp_and_function](std::string const& name,
                                             unsigned const num_comp,
                                             auto&& double_vec_from_loc_asm)
        {
            EXPECT_FALSE(map_name_to_num_comp_and_function.contains(name));
            map_name_to_num_comp_and_function[name] = {
                num_comp, std::move(double_vec_from_loc_asm)};
        });

    // checks //////////////////////////////////////////////////////////////////

    auto check = [&map_name_to_num_comp_and_function, &loc_asm](
                     std::string const& name,
                     unsigned const num_comp_expected,
                     std::vector<double> const& values_expected)
    {
        auto const it = map_name_to_num_comp_and_function.find(name);

        ASSERT_NE(map_name_to_num_comp_and_function.end(), it)
            << "No accessor found for ip data with name '" << name << "'";

        auto const& [num_comp, fct] = it->second;

        EXPECT_EQ(num_comp_expected, num_comp)
            << "Number of components differs for ip data with name '" << name
            << "'";
        EXPECT_THAT(fct(loc_asm),
                    testing::Pointwise(testing::DoubleEq(), values_expected))
            << "Values differ for ip data with name '" << name << "'";
    };

    // level 0
    check("scalar", 1, ref.scalar);
    check("vector", dim, ref.vector);
    check("kelvin", kv_size, ref.kelvin);

    // level 1
    check("scalar1", 1, ref.scalar1);
    check("vector1", dim, ref.vector1);
    check("kelvin1", kv_size, ref.kelvin1);

    // level 3
    check("scalar3", 1, ref.scalar3);
    check("vector3", dim, ref.vector3);
    check("kelvin3", kv_size, ref.kelvin3);
    check("matrix3", dim * 4, ref.matrix3);
    check("matrix3_1", dim * 4, ref.matrix3_1);
    check("matrix3_2", 4, ref.matrix3_2);

    // b levels
    check("scalar1b", 1, ref.scalar1b);
    check("scalar2b", 1, ref.scalar2b);
    check("scalar3b", 1, ref.scalar3b);
}

TYPED_TEST(ProcessLib_ReflectIPData, WriteTest)
{
    constexpr int dim = TypeParam::value;

    using LocAsm = LocAsmIF<dim>;

    std::size_t const num_int_pts = 8;
    LocAsm loc_asm(num_int_pts);

    auto const ref = ReferenceData<dim>::create(loc_asm, false);

    // data getters - used for checks //////////////////////////////////////////

    std::map<std::string, NumCompAndFunction<dim>>
        map_name_to_num_comp_and_function;

    ProcessLib::Reflection::forEachReflectedFlattenedIPDataAccessor<dim,
                                                                    LocAsm>(
        LocAsm::reflect(),
        [&map_name_to_num_comp_and_function](std::string const& name,
                                             unsigned const num_comp,
                                             auto&& double_vec_from_loc_asm)
        {
            EXPECT_FALSE(map_name_to_num_comp_and_function.contains(name));
            map_name_to_num_comp_and_function[name] = {
                num_comp, std::move(double_vec_from_loc_asm)};
        });

    // checks //////////////////////////////////////////////////////////////////

    auto check = [&map_name_to_num_comp_and_function, &loc_asm](
                     std::string const& name, std::size_t const size_expected,
                     std::vector<double> const& values_plain)
    {
        auto const it = map_name_to_num_comp_and_function.find(name);

        ASSERT_NE(map_name_to_num_comp_and_function.end(), it)
            << "No accessor found for ip data with name '" << name << "'";

        auto const& [num_comp, fct] = it->second;

        EXPECT_THAT(fct(loc_asm), testing::Each(testing::IsNan()))
            << "All values must be initialize to NaN in this unit test. Check "
               "failed for ip data with name '"
            << name << "'";

        auto const size = ProcessLib::Reflection::reflectSetIPData<dim>(
            name, values_plain.data(), loc_asm.ip_data_level1);

        EXPECT_EQ(size_expected, size)
            << "Unexpected size obtained for ip data with name '" << name
            << "'";

        // check set values via round-trip with getter tested in previous unit
        // test
        EXPECT_THAT(fct(loc_asm),
                    testing::Pointwise(testing::DoubleEq(), values_plain))
            << "Values not set correctly for ip data with name '" << name
            << "'";
    };

    check("scalar1", num_int_pts, ref.scalar1);
    check("vector1", num_int_pts, ref.vector1);
    check("kelvin1", num_int_pts, ref.kelvin1);

    check("scalar3", num_int_pts, ref.scalar3);
    check("vector3", num_int_pts, ref.vector3);
    check("kelvin3", num_int_pts, ref.kelvin3);
    check("matrix3", num_int_pts, ref.matrix3);
    check("matrix3_1", num_int_pts, ref.matrix3_1);
    check("matrix3_2", num_int_pts, ref.matrix3_2);

    check("scalar2b", num_int_pts, ref.scalar2b);
    check("scalar3b", num_int_pts, ref.scalar3b);
}

TYPED_TEST(ProcessLib_ReflectIPData, RawDataTypes)
{
    constexpr int dim = TypeParam::value;
    constexpr int kv_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(dim);

    namespace PRD = ProcessLib::Reflection::detail;

    // scalars
    static_assert(PRD::is_raw_data<double>::value);

    // vectors
    static_assert(PRD::is_raw_data<Eigen::Vector<double, dim>>::value);
    static_assert(PRD::is_raw_data<Eigen::RowVector<double, dim>>::value);

    // Kelvin vectors
    static_assert(PRD::is_raw_data<Eigen::Vector<double, kv_size>>::value);

    // matrices
    static_assert(PRD::is_raw_data<
                  Eigen::Matrix<double, dim, dim, Eigen::RowMajor>>::value);
    static_assert(PRD::is_raw_data<
                  Eigen::Matrix<double, dim, kv_size, Eigen::RowMajor>>::value);
    static_assert(PRD::is_raw_data<
                  Eigen::Matrix<double, kv_size, dim, Eigen::RowMajor>>::value);
    static_assert(
        PRD::is_raw_data<
            Eigen::Matrix<double, kv_size, kv_size, Eigen::RowMajor>>::value);

    // column major matrices are not supported in order to avoid confusion of
    // storage order
    static_assert(!PRD::is_raw_data<Eigen::Matrix<double, dim, dim>>::value);
    static_assert(
        !PRD::is_raw_data<Eigen::Matrix<double, dim, kv_size>>::value);
    static_assert(
        !PRD::is_raw_data<Eigen::Matrix<double, kv_size, dim>>::value);
    static_assert(
        !PRD::is_raw_data<Eigen::Matrix<double, kv_size, kv_size>>::value);
}
