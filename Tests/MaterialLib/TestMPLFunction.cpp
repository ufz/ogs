// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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

    static constexpr std::array<double, 3> coords = {10, 20, 35};
    ParameterLib::SpatialPosition const pos =
        ParameterLib::SpatialPosition{{0}, {0}, MathLib::Point3d{coords}};
    static constexpr double t = 10.;
};

TEST_F(MPLFunction, ScalarScalar)
{
    vars.temperature = 2.;

    MPL::Property const& f = MPL::Function{
        "test_function", {"temperature"}, {{"temperature", {"1"}}}, {}};

    ASSERT_EQ(2., f.value<double>(vars, pos, nan, nan));
    ASSERT_EQ(
        1., f.dValue<double>(vars, MPL::Variable::temperature, pos, nan, nan));

    MPL::Property const& f_t =
        MPL::Function{"function_t", {"t"}, {{"temperature", {"0"}}}, {}};
    ASSERT_EQ(10., f_t.value<double>(vars, pos, t, nan));
    ASSERT_EQ(
        0., f_t.dValue<double>(vars, MPL::Variable::temperature, pos, t, nan));

    MPL::Property const& f_temperature_t =
        MPL::Function{"function_temperature*t",
                      {"temperature * t"},
                      {{"temperature", {"t"}}},
                      {}};
    ASSERT_EQ(20., f_temperature_t.value<double>(vars, pos, t, nan));
    ASSERT_EQ(t,
              f_temperature_t.dValue<double>(vars, MPL::Variable::temperature,
                                             pos, t, nan));

    MPL::Property const& f_pos = MPL::Function{
        "function_x+y+z", {"x+y+z"}, {{"temperature", {"0"}}}, {}};
    ASSERT_EQ(coords[0] + coords[1] + coords[2],
              f_pos.value<double>(vars, pos, t, nan));
    ASSERT_EQ(
        0, f_pos.dValue<double>(vars, MPL::Variable::temperature, pos, t, nan));

    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;
    curves["linear_curve"] =
        std::make_unique<MathLib::PiecewiseLinearInterpolation>(
            std::vector<double>{0, 4}, std::vector<double>{0, 1}, true);
    MPL::Property const& f_curve =
        MPL::Function{"function_curve",
                      {"linear_curve(temperature)"},
                      {{"temperature", {"linear_curve(t)"}}},
                      curves};
    ASSERT_EQ(0.5, f_curve.value<double>(vars, pos, t, nan));
    ASSERT_EQ(
        1.,
        f_curve.dValue<double>(vars, MPL::Variable::temperature, pos, t, nan));
}

TEST_F(MPLFunction, ScalarVector)
{
    vars.temperature = 2.;

    MPL::Property const& f =
        MPL::Function{"test_function",
                      {"temperature", "temperature^2"},
                      {{"temperature", {"1", "2*temperature"}}},
                      {}};
    ASSERT_EQ((Eigen::Vector2d{2., 4.}),
              (f.value<Eigen::Vector2d>(vars, pos, t, nan)));
    ASSERT_EQ((Eigen::Vector2d{1., 4.}),
              f.dValue<Eigen::Vector2d>(vars, MPL::Variable::temperature, pos,
                                        t, nan));

    MPL::Property const& f_t =
        MPL::Function{"test_function",
                      {"temperature + t", "temperature^2-t"},
                      {{"temperature", {"1+t", "t^2*temperature"}}},
                      {}};
    ASSERT_EQ((Eigen::Vector2d{12., -6.}),
              (f_t.value<Eigen::Vector2d>(vars, pos, t, nan)));
    ASSERT_EQ((Eigen::Vector2d{11., 200.}),
              f_t.dValue<Eigen::Vector2d>(vars, MPL::Variable::temperature, pos,
                                          t, nan));

    MPL::Property const& f_pos =
        MPL::Function{"test_function",
                      {"temperature + y", "temperature^2-z"},
                      {{"temperature", {"1+x+y-z", "y^2*temperature"}}},
                      {}};
    ASSERT_EQ((Eigen::Vector2d{22., -31.}),
              (f_pos.value<Eigen::Vector2d>(vars, pos, t, nan)));
    ASSERT_EQ((Eigen::Vector2d{-4, 800.}),
              f_pos.dValue<Eigen::Vector2d>(vars, MPL::Variable::temperature,
                                            pos, t, nan));

    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;
    curves["linear_curve"] =
        std::make_unique<MathLib::PiecewiseLinearInterpolation>(
            std::vector<double>{0, 4}, std::vector<double>{0, 1}, true);
    MPL::Property const& f_curve = MPL::Function{
        "function_curve",
        {"linear_curve(temperature)", "linear_curve(temperature)^2"},
        {{"temperature", {"linear_curve(t)", "0"}}},
        curves};
    ASSERT_EQ((Eigen::Vector2d{0.5, 0.25}),
              (f_curve.value<Eigen::Vector2d>(vars, pos, t, nan)));
    ASSERT_EQ((Eigen::Vector2d{1, 0}),
              f_curve.dValue<Eigen::Vector2d>(vars, MPL::Variable::temperature,
                                              pos, t, nan));
}

// Regression test pinning the row-major layout of matrix-valued Function
// output. A 4-expression value maps to Eigen::Matrix2d in row major order (see
// fromArray in MaterialLib/MPL/Property.h, consistent with fromVector).
// Distinct asymmetric literals make a column-major regression observable:
// column major would yield [[1,3],[2,4]] (m(0,1) == 3 instead of 2).
TEST_F(MPLFunction, ScalarMatrix2x2)
{
    MPL::Property const& f =
        MPL::Function{"test_function", {"1", "2", "3", "4"}, {}, {}};

    Eigen::Matrix2d const expected =
        (Eigen::Matrix2d() << 1, 2, 3, 4).finished();
    ASSERT_EQ(expected, (f.value<Eigen::Matrix2d>(vars, pos, t, nan)));
}

// As ScalarMatrix2x2, but for the separate 9-element branch mapping to
// Eigen::Matrix3d. Column-major regression would give m(0,1) == 4 instead of 2.
TEST_F(MPLFunction, ScalarMatrix3x3)
{
    MPL::Property const& f = MPL::Function{
        "test_function", {"1", "2", "3", "4", "5", "6", "7", "8", "9"}, {}, {}};

    Eigen::Matrix3d const expected =
        (Eigen::Matrix3d() << 1, 2, 3, 4, 5, 6, 7, 8, 9).finished();
    ASSERT_EQ(expected, (f.value<Eigen::Matrix3d>(vars, pos, t, nan)));
}

TEST_F(MPLFunction, KelvinVector2Scalar)
{
    vars.stress.emplace<KV2>(1, 2, 3, 4 * std::sqrt(2.));
    vars.temperature = 273.15;

    MPL::Property const& f =
        MPL::Function{"test_function", {"avg(stress) * temperature"}, {}, {}};
    ASSERT_EQ(((1 + 4) * 0.5) * 273.15, f.value<double>(vars, pos, t, nan));

    MPL::Property const& f_t = MPL::Function{
        "test_function", {"avg(stress) * temperature + t"}, {}, {}};
    ASSERT_EQ(((1 + 4) * 0.5) * 273.15 + 10,
              f_t.value<double>(vars, pos, t, nan));

    MPL::Property const& f_pos = MPL::Function{
        "test_function", {"avg(stress) * temperature + x + y + z"}, {}, {}};
    ASSERT_EQ(((1 + 4) * 0.5) * 273.15 + coords[0] + coords[1] + coords[2],
              f_pos.value<double>(vars, pos, t, nan));

    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;
    curves["linear_curve"] =
        std::make_unique<MathLib::PiecewiseLinearInterpolation>(
            std::vector<double>{0, 5}, std::vector<double>{10, 0}, true);

    MPL::Property const& f_curve = MPL::Function{
        "test_function", {"linear_curve(avg(stress))"}, {}, curves};
    ASSERT_EQ(5, f_curve.value<double>(vars, pos, t, nan));
}

TEST_F(MPLFunction, KelvinVector3Scalar)
{
    vars.stress.emplace<KV3>(1, 2, 3, 4 * std::sqrt(2.), 5 * std::sqrt(2.),
                             6 * std::sqrt(2.));
    vars.temperature = 273.15;

    MPL::Property const& f =
        MPL::Function{"test_function", {"avg(stress) * temperature"}, {}, {}};
    ASSERT_EQ(((1 + 6) * 0.5) * 273.15, f.value<double>(vars, pos, t, nan));

    MPL::Property const& f_t = MPL::Function{
        "test_function", {"avg(stress) * temperature + t"}, {}, {}};
    ASSERT_EQ(((1 + 6) * 0.5) * 273.15 + 10,
              f_t.value<double>(vars, pos, t, nan));

    MPL::Property const& f_pos = MPL::Function{
        "test_function", {"avg(stress) * temperature + x + y + z"}, {}, {}};
    ASSERT_EQ(((1 + 6) * 0.5) * 273.15 + coords[0] + coords[1] + coords[2],
              f_pos.value<double>(vars, pos, t, nan));

    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;
    curves["linear_curve"] =
        std::make_unique<MathLib::PiecewiseLinearInterpolation>(
            std::vector<double>{0, 100}, std::vector<double>{0, 10}, true);

    MPL::Property const& f_curve = MPL::Function{
        "test_function", {"linear_curve(avg(stress) + x)"}, {}, curves};
    ASSERT_EQ(1.35, f_curve.value<double>(vars, pos, t, nan));
}

TEST_F(MPLFunction, KelvinVector23Scalar)
{
    // Mixed dimensions of vectorial quantities are expected to fail.
    vars.stress.emplace<KV3>(1, 2, 3, 4, 5, 6);
    vars.total_strain.emplace<KV2>(1, 2, 3, 4);

    MPL::Property const& f = MPL::Function{
        "test_function", {"avg(stress) * avg(total_strain)"}, {}, {}};
    ASSERT_ANY_THROW(f.value<double>(vars, pos, t, nan));

    MPL::Property const& f_t = MPL::Function{
        "test_function", {"avg(stress) * avg(total_strain) - t"}, {}, {}};
    ASSERT_ANY_THROW(f_t.value<double>(vars, pos, t, nan));

    MPL::Property const& f_pos =
        MPL::Function{"test_function",
                      {"avg(stress) * avg(total_strain) + x + y + z"},
                      {},
                      {}};
    ASSERT_ANY_THROW(f_pos.value<double>(vars, pos, t, nan));

    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;
    curves["linear_curve"] =
        std::make_unique<MathLib::PiecewiseLinearInterpolation>(
            std::vector<double>{0, 100}, std::vector<double>{0, 10}, true);
    MPL::Property const& f_curve =
        MPL::Function{"test_function",
                      {"avg(stress) * avg(total_strain) + linear_curve(x)"},
                      {},
                      curves};
    ASSERT_ANY_THROW(f_curve.value<double>(vars, pos, t, nan));
}

TEST_F(MPLFunction, VectorizedTensor2Scalar)
{
    vars.deformation_gradient.emplace<VT2>(1, 2, 3, 4, 5);
    vars.temperature = 273.15;

    MPL::Property const& f = MPL::Function{
        "test_function", {"avg(deformation_gradient) * temperature"}, {}, {}};
    ASSERT_EQ(((1 + 5) * 0.5) * 273.15, f.value<double>(vars, pos, t, nan));

    MPL::Property const& f_t =
        MPL::Function{"test_function",
                      {"t * avg(deformation_gradient) * temperature"},
                      {},
                      {}};
    ASSERT_EQ(10 * ((1 + 5) * 0.5) * 273.15,
              f_t.value<double>(vars, pos, t, nan));

    MPL::Property const& f_pos =
        MPL::Function{"test_function",
                      {"avg(deformation_gradient) * temperature + x+ y + z"},
                      {},
                      {}};
    ASSERT_EQ(((1 + 5) * 0.5) * 273.15 + coords[0] + coords[1] + coords[2],
              f_pos.value<double>(vars, pos, t, nan));

    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;
    curves["linear_curve"] =
        std::make_unique<MathLib::PiecewiseLinearInterpolation>(
            std::vector<double>{0, 100}, std::vector<double>{0, 10}, true);
    MPL::Property const& f_curve =
        MPL::Function{"test_function",
                      {"avg(deformation_gradient) * linear_curve(y)"},
                      {},
                      curves};
    ASSERT_EQ(6, f_curve.value<double>(vars, pos, t, nan));
}

TEST_F(MPLFunction, VectorizedTensor3Scalar)
{
    vars.deformation_gradient.emplace<VT3>(1, 2, 3, 4, 5, 6, 7, 8, 9);
    vars.temperature = 273.15;

    MPL::Property const& f = MPL::Function{
        "test_function", {"avg(deformation_gradient) * temperature"}, {}, {}};
    ASSERT_EQ(((1 + 9) * 0.5) * 273.15, f.value<double>(vars, pos, nan, nan));

    MPL::Property const& f_pos =
        MPL::Function{"test_function",
                      {"avg(deformation_gradient) * temperature/x + y + z"},
                      {},
                      {}};
    ASSERT_EQ(((1 + 9) * 0.5) * 273.15 / coords[0] + coords[1] + coords[2],
              f_pos.value<double>(vars, pos, nan, nan));

    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;
    curves["linear_curve"] =
        std::make_unique<MathLib::PiecewiseLinearInterpolation>(
            std::vector<double>{0, 100}, std::vector<double>{0, 10}, true);
    MPL::Property const& f_curve = MPL::Function{
        "test_function",
        {"avg(deformation_gradient) * linear_curve(temperature/x)"},
        {},
        curves};
    ASSERT_EQ(13.6575, f_curve.value<double>(vars, pos, nan, nan));
}

TEST_F(MPLFunction, ScalarUninitialized)
{
    // The vars.temperature = is not initialized.

    MPL::Property const& f = MPL::Function{
        "test_function", {"temperature"}, {{"temperature", {"1"}}}, {}};
    ASSERT_ANY_THROW(f.value<double>(vars, {}, nan, nan));
    ASSERT_ANY_THROW(
        f.dValue<double>(vars, MPL::Variable::temperature, {}, nan, nan));
}

TEST_F(MPLFunction, KelvinVector2Uninitialized)
{
    vars.temperature = 273.15;
    // The vars.stress is not initialized.

    MPL::Property const& f =
        MPL::Function{"test_function", {"avg(stress) * temperature"}, {}, {}};
    ASSERT_ANY_THROW(f.value<double>(vars, {}, nan, nan));
}

TEST_F(MPLFunction, VectorizedTensor2Uninitialized)
{
    vars.temperature = 273.15;
    // The vars.deformation_gradient is not initialized.

    MPL::Property const& f = MPL::Function{
        "test_function", {"avg(deformation_gradient) * temperature"}, {}, {}};
    ASSERT_ANY_THROW(f.value<double>(vars, {}, nan, nan));
}

TEST_F(MPLFunction, KelvinVector3Uninitialized)
{
    vars.temperature = 273.15;
    // The vars.stress is not initialized.

    MPL::Property const& f =
        MPL::Function{"test_function", {"avg(stress) * temperature"}, {}, {}};
    ASSERT_ANY_THROW(f.value<double>(vars, {}, nan, nan));
}

TEST_F(MPLFunction, VectorizedTensor3Uninitialized)
{
    vars.temperature = 273.15;
    // The vars.deformation_gradient is not initialized.

    MPL::Property const& f = MPL::Function{
        "test_function", {"avg(deformation_gradient) * temperature"}, {}, {}};
    ASSERT_ANY_THROW(f.value<double>(vars, {}, nan, nan));
}

// ============================================================================
// OpenMP Thread-Safety Tests
// ============================================================================

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Tests/AutoCheckTools.h"

namespace ac = autocheck;

struct TestInput
{
    double t;
    MathLib::Point3d pos;
    MathLib::KelvinVector::KelvinVectorType<3> stress;
    double temperature;
};

std::ostream& operator<<(std::ostream& os, TestInput const& input)
{
    os << "TestInput{t=" << input.t << ", pos=[" << input.pos[0] << ","
       << input.pos[1] << "," << input.pos[2] << "], stress=["
       << input.stress[0] << "," << input.stress[1] << "," << input.stress[2]
       << "," << input.stress[3] << "," << input.stress[4] << ","
       << input.stress[5] << "], temperature=" << input.temperature << "}";
    return os;
}

struct GeneratorForTestInput
{
    using result_type = TestInput;

    result_type operator()(std::size_t size)
    {
        ac::generator<double> gen;
        MathLib::KelvinVector::KelvinVectorType<3> stress;
        for (int j = 0; j < 6; ++j)
        {
            stress[j] = gen(size);
        }
        return {
            gen(size) * 2.0,
            MathLib::Point3d{{gen(size), gen(size), gen(size)}},
            stress,
            gen(size),
        };
    }
};

struct GeneratorForTestInputVector
{
    using result_type = std::vector<TestInput>;

    result_type operator()(std::size_t size)
    {
        std::size_t const n = size > 0 ? size : 20;
        result_type result;
        result.reserve(n);
        ac::generator<double> gen;
        for (std::size_t i = 0; i < n; ++i)
        {
            MathLib::KelvinVector::KelvinVectorType<3> stress;
            for (int j = 0; j < 6; ++j)
            {
                stress[j] = gen(size);
            }
            // Unlike GeneratorForTestInput, t is scaled by 1/n so that the
            // generated times stay within the sloped [0, 2] range of the test
            // curve (rather than saturating it) as n grows.
            result.push_back({
                gen(size) * 2.0 / static_cast<double>(n),
                MathLib::Point3d{{gen(size), gen(size), gen(size)}},
                stress,
                gen(size),
            });
        }
        return result;
    }
};

struct MPLFunctionOpenMPTest : public ::testing::Test,
                               public ::testing::WithParamInterface<int>
{
    static double constexpr nan = std::numeric_limits<double>::quiet_NaN();

    void SetUp() override
    {
#ifdef _OPENMP
        original_max_threads = omp_get_max_threads();
        omp_set_num_threads(GetParam());
#endif
        curves["linear"] = createLinearTestCurve();
    }

    void TearDown() override
    {
#ifdef _OPENMP
        omp_set_num_threads(original_max_threads);
#endif
    }

    static std::unique_ptr<MathLib::PiecewiseLinearInterpolation>
    createLinearTestCurve()
    {
        std::vector<double> support_points{0.0, 1.0, 2.0};
        std::vector<double> values{0.0, 1.0, 2.0};
        return std::make_unique<MathLib::PiecewiseLinearInterpolation>(
            std::move(support_points), std::move(values));
    }

    int original_max_threads = 1;
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;
    ac::gtest_reporter gtest_reporter;
};

TEST_P(MPLFunctionOpenMPTest, DeterministicSameInputEvaluation)
{
    auto property = [&](TestInput const& input)
    {
    // Thread-safety check: a single Function instance is evaluated
    // concurrently on every thread with one and the same input. If the
    // evaluation has no data races (no shared mutable state), all threads
    // must agree on the result, and that common result must match the
    // analytically known value. This is verified separately for value()
    // and dValue().
#ifdef _OPENMP
        int const num_threads = omp_get_max_threads();
#else
        int const num_threads = 1;
#endif

        MPL::Property const& f = MPL::Function{"test",
                                               {"temperature + t + x + y + z"},
                                               {{"temperature", {"1"}}},
                                               {}};

        MPL::VariableArray vars;
        vars.temperature = input.temperature;

        ParameterLib::SpatialPosition pos;
        pos.setCoordinates(input.pos);

        std::vector<double> value_results(num_threads);

#pragma omp parallel for
        for (int i = 0; i < num_threads; ++i)
        {
            // all testees are called with the same inputs...
            value_results[i] = f.value<double>(vars, pos, input.t, nan);
        }

        for (int i = 1; i < num_threads; ++i)
        {
            // ... and must return the same results
            if (value_results[i] != value_results[0])
            {
                return false;
            }
        }

        std::vector<double> dvalue_results(num_threads);

#pragma omp parallel for
        for (int i = 0; i < num_threads; ++i)
        {
            // same for the derivative: identical inputs on every thread...
            dvalue_results[i] = f.dValue<double>(
                vars, MPL::Variable::temperature, pos, input.t, nan);
        }

        for (int i = 1; i < num_threads; ++i)
        {
            // ... must yield identical results
            if (dvalue_results[i] != dvalue_results[0])
            {
                return false;
            }
        }

        // Threads agree (checked above); now confirm correctness against the
        // closed-form expression "temperature + t + x + y + z" and its
        // temperature-derivative (constant 1.0).
        double const expected_value = input.temperature + input.t +
                                      input.pos[0] + input.pos[1] +
                                      input.pos[2];
        double const expected_dvalue = 1.0;

        return std::abs(value_results[0] - expected_value) < 1e-10 &&
               std::abs(dvalue_results[0] - expected_dvalue) < 1e-10;
    };

    ac::check<TestInput>(property, 100,
                         ac::make_arbitrary(GeneratorForTestInput{}),
                         gtest_reporter);
}

TEST_P(MPLFunctionOpenMPTest, SequentialParallelEquivalence)
{
    auto property = [&](TestInput const& input_base)
    {
        std::vector<TestInput> inputs;
        inputs.reserve(20);
        for (int i = 0; i < 20; ++i)
        {
            inputs.push_back({
                input_base.t + i * 0.05,
                input_base.pos,
                input_base.stress,
                input_base.temperature + i * 0.1,
            });
        }

        MPL::Property const& f_seq = MPL::Function{
            "test", {"temperature * t + x * y"}, {{"temperature", {"t"}}}, {}};
        MPL::Property const& f_par = MPL::Function{
            "test", {"temperature * t + x * y"}, {{"temperature", {"t"}}}, {}};

        std::vector<std::pair<double, double>> seq_results;
        seq_results.reserve(inputs.size());
        for (auto const& inp : inputs)
        {
            MPL::VariableArray vars;
            vars.temperature = inp.temperature;

            ParameterLib::SpatialPosition pos;
            pos.setCoordinates(inp.pos);

            double val = f_seq.value<double>(vars, pos, inp.t, nan);
            double dval = f_seq.dValue<double>(vars, MPL::Variable::temperature,
                                               pos, inp.t, nan);
            seq_results.emplace_back(val, dval);
        }

        std::vector<std::pair<double, double>> par_results(inputs.size());

        int const num_threads = static_cast<int>(inputs.size());
#pragma omp parallel for
        for (int i = 0; i < num_threads; ++i)
        {
            auto const& inp = inputs[i];
            MPL::VariableArray vars;
            vars.temperature = inp.temperature;

            ParameterLib::SpatialPosition pos;
            pos.setCoordinates(inp.pos);

            double val = f_par.value<double>(vars, pos, inp.t, nan);
            double dval = f_par.dValue<double>(vars, MPL::Variable::temperature,
                                               pos, inp.t, nan);
            par_results[i] = {val, dval};
        }

        return seq_results == par_results;
    };

    ac::check<TestInput>(property, 100,
                         ac::make_arbitrary(GeneratorForTestInput{}),
                         gtest_reporter);
}

TEST_P(MPLFunctionOpenMPTest, IndependentThreadIsolationWithCurves)
{
    auto property = [&](std::vector<TestInput> const& inputs)
    {
        MPL::Property const& f = MPL::Function{"test",
                                               {"linear(t) + temperature + x"},
                                               {{"temperature", {"linear(t)"}}},
                                               curves};

        std::vector<double> value_results(inputs.size());

        int const num_threads = static_cast<int>(inputs.size());
#pragma omp parallel for
        for (int i = 0; i < num_threads; ++i)
        {
            auto const& inp = inputs[i];
            MPL::VariableArray vars;
            vars.temperature = inp.temperature;

            ParameterLib::SpatialPosition pos;
            pos.setCoordinates(inp.pos);

            value_results[i] = f.value<double>(vars, pos, inp.t, nan);
        }

        std::vector<double> dvalue_results(inputs.size());

#pragma omp parallel for
        for (int i = 0; i < num_threads; ++i)
        {
            auto const& inp = inputs[i];
            MPL::VariableArray vars;
            vars.temperature = inp.temperature;

            ParameterLib::SpatialPosition pos;
            pos.setCoordinates(inp.pos);

            dvalue_results[i] = f.dValue<double>(
                vars, MPL::Variable::temperature, pos, inp.t, nan);
        }

        auto linear = [](double t) -> double
        {
            if (t < 0.0)
            {
                return 0.0;
            }
            if (t > 2.0)
            {
                return 2.0;
            }
            return t;
        };

        for (std::size_t i = 0; i < inputs.size(); ++i)
        {
            auto const& inp = inputs[i];
            double const expected_value =
                linear(inp.t) + inp.temperature + inp.pos[0];
            if (std::abs(value_results[i] - expected_value) > 1e-10)
            {
                return false;
            }

            double const expected_dvalue = linear(inp.t);
            if (std::abs(dvalue_results[i] - expected_dvalue) > 1e-10)
            {
                return false;
            }
        }

        return true;
    };

    ac::check<std::vector<TestInput>>(
        property, 100, ac::make_arbitrary(GeneratorForTestInputVector{}),
        gtest_reporter);
}

#ifdef _OPENMP
INSTANTIATE_TEST_SUITE_P(MPLFunctionOpenMP,
                         MPLFunctionOpenMPTest,
                         ::testing::Values(1, 2, 4, 8, 16, 32, 64),
                         [](testing::TestParamInfo<int> const& info)
                         { return "threads_" + std::to_string(info.param); });
#else
INSTANTIATE_TEST_SUITE_P(MPLFunctionOpenMP,
                         MPLFunctionOpenMPTest,
                         ::testing::Values(1),
                         [](testing::TestParamInfo<int> const& info)
                         { return "threads_" + std::to_string(info.param); });
#endif
