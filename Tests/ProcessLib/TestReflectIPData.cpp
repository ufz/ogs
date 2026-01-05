// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <numeric>
#include <random>

#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"
#include "ProcessLib/Output/CellAverageAlgorithm.h"
#include "ProcessLib/Reflection/ReflectionForIPWriters.h"
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
using Kelvin1Data =
    BaseLib::StrongType<MathLib::KelvinVector::KelvinVectorType<Dim>,
                        struct Kelvin1DataTag>;

constexpr std::string_view ioName(struct Kelvin1DataTag*)
{
    return "kelvin1";
}

template <int Dim>
struct Level1
{
    Kelvin1Data<Dim> kelvin1;
    Eigen::Vector<double, Dim> vector1;
    double scalar1;
    Level2<Dim> level2;
    Level2b<Dim> level2b;

    static auto reflect()
    {
        using namespace ProcessLib::Reflection;
        return std::tuple{makeReflectionData(&Level1::kelvin1),
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

using ScalarData = BaseLib::StrongType<double, struct ScalarDataTag>;

constexpr std::string_view ioName(struct ScalarDataTag*)
{
    return "scalar";
}

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

    std::vector<ScalarData> ip_data_scalar;
    std::vector<Eigen::Vector<double, Dim>> ip_data_vector;
    std::vector<MathLib::KelvinVector::KelvinVectorType<Dim>> ip_data_kelvin;
    std::vector<Level1<Dim>> ip_data_level1;
    std::vector<Level1b<Dim>> ip_data_level1b;

    static auto reflect()
    {
        using namespace ProcessLib::Reflection;
        return std::tuple{
            makeReflectionData(&LocAsmIF::ip_data_scalar),
            makeReflectionData("vector", &LocAsmIF::ip_data_vector),
            makeReflectionData("kelvin", &LocAsmIF::ip_data_kelvin),
            makeReflectionData(&LocAsmIF::ip_data_level1),
            makeReflectionData(&LocAsmIF::ip_data_level1b)};
    }

    static auto getReflectionDataForOutput() { return reflect(); }
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
                // the local assembler stores Kelvin vector data...
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
    // ... the reference data used in the tests is not Kelvin-mapped!
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
private:
    ReferenceData() = default;

public:
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

    // Computes reference data and initializes the internal (integration point)
    // data of the passed \c loc_asm.
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
                return *loc_asm.ip_data_scalar[ip];
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
                return *loc_asm.ip_data_level1[ip].kelvin1;
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

template <int dim>
void checkAll(ReferenceData<dim> const& ref, auto&& checker,
              bool const check_level_0_data = true)
{
    auto constexpr kv_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(dim);

    // The storage of "level 0 data" (i.e., non-reflected data directly in the
    // local assembler, e.g. a std::vector<double> of integration point data of
    // a scalar value) is not compatible with the implemented IP data input
    // logic. Therefore we don't test level 0 data in that case.
    // This incompatibility is not a problem in practice, as such data is not
    // used ATM in OGS's local assemblers.
    if (check_level_0_data)
    {
        checker("scalar", 1, ref.scalar);
        checker("vector", dim, ref.vector);
        checker("kelvin", kv_size, ref.kelvin);
    }

    // level 1
    checker("scalar1", 1, ref.scalar1);
    checker("vector1", dim, ref.vector1);
    checker("kelvin1", kv_size, ref.kelvin1);

    // level 3
    checker("scalar3", 1, ref.scalar3);
    checker("vector3", dim, ref.vector3);
    checker("kelvin3", kv_size, ref.kelvin3);
    checker("matrix3", dim * 4, ref.matrix3);
    checker("matrix3_1", dim * 4, ref.matrix3_1);
    checker("matrix3_2", 4, ref.matrix3_2);

    // b levels
    checker("scalar1b", 1, ref.scalar1b);
    checker("scalar2b", 1, ref.scalar2b);
    checker("scalar3b", 1, ref.scalar3b);
}

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

    checkAll(ref, check);
}

TYPED_TEST(ProcessLib_ReflectIPData, WriteTest)
{
    constexpr int dim = TypeParam::value;

    using LocAsm = LocAsmIF<dim>;

    std::size_t const num_int_pts = 8;
    LocAsm loc_asm(num_int_pts);

    auto const ref = ReferenceData<dim>::create(loc_asm, false);

    // set up data getters - used for checks ///////////////////////////////////
    // that critically relies on the read test above being successful!

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
                     std::string const& name, unsigned const num_comp,
                     std::vector<double> const& values_plain)
    {
        auto const it = map_name_to_num_comp_and_function.find(name);

        ASSERT_NE(map_name_to_num_comp_and_function.end(), it)
            << "No accessor found for ip data with name '" << name << "'";

        auto const& [num_comp_2, fct] = it->second;

        // consistency checks
        ASSERT_EQ(num_comp, num_comp_2);
        ASSERT_EQ(0, values_plain.size() % num_comp);
        auto const num_int_pts_expected = values_plain.size() / num_comp;

        EXPECT_THAT(fct(loc_asm), testing::Each(testing::IsNan()))
            << "All values must be initialized to NaN in this unit test. Check "
               "failed for ip data with name '"
            << name << "'";

        // function under test /////////////////////////////////////////////////

        auto const num_int_pts_actual =
            ProcessLib::Reflection::reflectSetIPData<dim>(
                name, values_plain.data(), loc_asm.ip_data_level1);
        auto const num_int_pts_actual_1 =
            ProcessLib::Reflection::reflectSetIPData<dim>(
                name, values_plain.data(), loc_asm.ip_data_level1b);

        // end function under test /////////////////////////////////////////////

        ASSERT_EQ(num_int_pts_expected, num_int_pts_actual)
            << "Unexpected number of integration points obtained for ip data "
               "with name '"
            << name << "'";
        ASSERT_EQ(num_int_pts_expected, num_int_pts_actual_1);

        // check set values via round-trip with getter tested in previous unit
        // test
        EXPECT_THAT(fct(loc_asm),
                    testing::Pointwise(testing::DoubleEq(), values_plain))
            << "Values not set correctly for ip data with name '" << name
            << "'";
    };

    checkAll(ref, check, false);
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

TYPED_TEST(ProcessLib_ReflectIPData, IPWriterTest)
{
    constexpr int dim = TypeParam::value;

    using LocAsm = LocAsmIF<dim>;

    std::size_t const num_int_pts = 8;
    std::vector<std::unique_ptr<LocAsm>> loc_asms;
    // emplace...new is a workaround for an MSVC compiler warning:
    // C:\Program Files (x86)\Microsoft Visual
    // Studio\2019\Community\VC\Tools\MSVC\14.29.30133\include\memory(3382):
    // warning C4244: 'argument': conversion from 'const uintptr_t' to 'const
    // unsigned int', possible loss of data
    // C:/Users/gitlab/gitlab/_b/bg4d5s_d/0/ogs/ogs/Tests/ProcessLib/
    // TestReflectIPData.cpp(792):
    // note: see reference to function template instantiation
    // 'std::unique_ptr<LocAsm,std::default_delete<LocAsm>>
    // std::make_unique<LocAsm,const size_t&,0>(const size_t &)' being compiled
    loc_asms.emplace_back(new LocAsm{num_int_pts});

    std::unique_ptr<MeshLib::Mesh> mesh{
        MeshToolsLib::MeshGenerator::generateLineMesh(1.0, 1)};

    auto const ref = ReferenceData<dim>::create(*loc_asms.front(), true);

    // function under test /////////////////////////////////////////////////////

    std::vector<std::unique_ptr<MeshLib::IntegrationPointWriter>>
        integration_point_writers;
    constexpr int integration_order = 0xc0ffee;  // dummy value
    ProcessLib::Reflection::addReflectedIntegrationPointWriters<dim>(
        LocAsm::getReflectionDataForOutput(), integration_point_writers,
        integration_order, loc_asms);

    // this finally invokes the IP writers
    MeshLib::addIntegrationPointDataToMesh(*mesh, integration_point_writers);

    // checks //////////////////////////////////////////////////////////////////

    auto check = [&mesh](std::string const& name,
                         unsigned const num_comp,
                         std::vector<double> const& values_expected)
    {
        auto const& props = mesh->getProperties();

        ASSERT_TRUE(props.existsPropertyVector<double>(name + "_ip"))
            << "Property vector '" << name + "_ip"
            << "' does not exist in the mesh.";

        ASSERT_EQ(0, values_expected.size() % num_comp);
        auto const num_int_pts = values_expected.size() / num_comp;

        auto const& ip_data_actual =
            *props.getPropertyVector<double>(name + "_ip");

        ASSERT_EQ(MeshLib::MeshItemType::IntegrationPoint,
                  ip_data_actual.getMeshItemType());
        ASSERT_EQ(mesh->getNumberOfElements() * num_int_pts,
                  ip_data_actual.getNumberOfTuples());
        ASSERT_EQ(num_comp, ip_data_actual.getNumberOfGlobalComponents());

        EXPECT_THAT(ip_data_actual,
                    testing::Pointwise(testing::DoubleEq(), values_expected))
            << "Values differ for ip data with name '" << name << "'";
    };

    checkAll(ref, check);
}

TYPED_TEST(ProcessLib_ReflectIPData, CellAverageTest)
{
    constexpr int dim = TypeParam::value;

    using LocAsm = LocAsmIF<dim>;

    std::size_t const num_int_pts = 8;
    std::vector<std::unique_ptr<LocAsm>> loc_asms;
    // emplace...new is a workaround for an MSVC compiler warning:
    // C:\Program Files (x86)\Microsoft Visual
    // Studio\2019\Community\VC\Tools\MSVC\14.29.30133\include\memory(3382):
    // warning C4244: 'argument': conversion from 'const uintptr_t' to 'const
    // unsigned int', possible loss of data
    // C:/Users/gitlab/gitlab/_b/bg4d5s_d/0/ogs/ogs/Tests/ProcessLib/
    // TestReflectIPData.cpp(792):
    // note: see reference to function template instantiation
    // 'std::unique_ptr<LocAsm,std::default_delete<LocAsm>>
    // std::make_unique<LocAsm,const size_t&,0>(const size_t &)' being compiled
    loc_asms.emplace_back(new LocAsm{num_int_pts});

    std::unique_ptr<MeshLib::Mesh> mesh{
        MeshToolsLib::MeshGenerator::generateLineMesh(1.0, 1)};

    auto const ref = ReferenceData<dim>::create(*loc_asms.front(), true);

    // function under test /////////////////////////////////////////////////////

    ProcessLib::CellAverageData cell_average_data{*mesh};
    computeCellAverages<dim>(cell_average_data, loc_asms);

    // checks //////////////////////////////////////////////////////////////////

    auto check = [&mesh](std::string const& name,
                         unsigned const num_comp_expected,
                         std::vector<double> const& ip_data)
    {
        ASSERT_EQ(num_int_pts * num_comp_expected, ip_data.size());

        // TODO the implementation here in the unit test computing cell average
        // reference data might be too similar to the production code. In fact,
        // it's almost the same code.

        // each integration point corresponds to a column in the mapped
        // matrix, vector components/matrix entries are stored contiguously
        // in memory.
        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                       Eigen::ColMajor>>
            ip_data_mat{ip_data.data(), num_comp_expected, num_int_pts};

        Eigen::VectorXd const cell_avg_expected = ip_data_mat.rowwise().mean();

        auto const& props = mesh->getProperties();

        ASSERT_TRUE(props.existsPropertyVector<double>(name + "_avg"))
            << "Property vector '" << name + "_avg"
            << "' does not exist in the mesh.";

        auto const& cell_avg_actual =
            *props.getPropertyVector<double>(name + "_avg");

        ASSERT_EQ(MeshLib::MeshItemType::Cell,
                  cell_avg_actual.getMeshItemType());
        ASSERT_EQ(mesh->getNumberOfElements(),
                  cell_avg_actual.getNumberOfTuples());
        ASSERT_EQ(num_comp_expected,
                  cell_avg_actual.getNumberOfGlobalComponents());

        EXPECT_THAT(cell_avg_actual,
                    testing::Pointwise(testing::DoubleEq(), cell_avg_expected))
            << "Values differ for cell average data with name '" << name << "'";
    };

    checkAll(ref, check);
}
