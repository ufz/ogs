// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ParameterLib/FunctionParameter.h"
#include "ParameterLib/SpatialPosition.h"
#include "Tests/AutoCheckTools.h"

namespace ac = autocheck;

static std::unique_ptr<MathLib::PiecewiseLinearInterpolation>
createLinearTestCurveFunctionParameter()
{
    std::vector<double> support_points{0.0, 1.0, 2.0};
    std::vector<double> values{0.0, 1.0, 2.0};
    return std::make_unique<MathLib::PiecewiseLinearInterpolation>(
        std::move(support_points), std::move(values));
}

struct GeneratorForVectorArray4
{
    using result_type = std::vector<std::array<double, 4>>;

    result_type operator()(std::size_t size)
    {
        std::size_t const n = size > 0 ? size : 20;
        result_type result;
        result.reserve(n);
        ac::generator<double> gen;
        for (std::size_t i = 0; i < n; ++i)
        {
            result.emplace_back(std::array<double, 4>{gen(size), gen(size),
                                                      gen(size), gen(size)});
        }
        return result;
    }
};

struct FunctionParameterOpenMPTest : public ::testing::Test,
                                     public ::testing::WithParamInterface<int>
{
    void SetUp() override
    {
#ifdef _OPENMP
        original_max_threads_ = omp_get_max_threads();
        omp_set_num_threads(GetParam());
#endif
        curves_["linear"] = createLinearTestCurveFunctionParameter();
    }

    void TearDown() override
    {
#ifdef _OPENMP
        omp_set_num_threads(original_max_threads_);
#endif
    }

    int original_max_threads_ = 1;
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves_;
    ac::gtest_reporter gtest_reporter;
};

TEST_P(FunctionParameterOpenMPTest, DeterministicSameInputEvaluation)
{
    auto property = [&](std::array<double, 4> const& input)
    {
#ifdef _OPENMP
        int const num_threads = omp_get_max_threads();
#else
        int const num_threads = 1;
#endif

        double const t = input[0];
        double const x = input[1];
        double const y = input[2];
        double const z = input[3];

        ParameterLib::FunctionParameter<double> param("test", {"x + y + z + t"},
                                                      {});
        ParameterLib::SpatialPosition pos;
        pos.setCoordinates(MathLib::Point3d{{x, y, z}});

        std::vector<std::vector<double>> results(
            static_cast<std::size_t>(num_threads));

#pragma omp parallel for
        for (int i = 0; i < num_threads; ++i)
        {
            results[static_cast<std::size_t>(i)] = param(t, pos);
        }

        for (int i = 1; i < num_threads; ++i)
        {
            if (results[static_cast<std::size_t>(i)] != results[0])
            {
                return false;
            }
        }

        double const expected = t + x + y + z;
        return std::abs(results[0][0] - expected) < 1e-10;
    };

    ac::IntervalGenerator<double> gen(-5.0, 5.0);

    ac::check<std::array<double, 4>>(
        property,
        100,
        ac::make_arbitrary(
            ac::randomTupleGenerator<double, 4, decltype(gen)>{gen}),
        gtest_reporter);
}

TEST_P(FunctionParameterOpenMPTest, SequentialParallelEquivalence)
{
    auto property = [&](std::array<double, 4> const& input)
    {
        std::vector<std::array<double, 4>> inputs;
        inputs.reserve(20);
        for (int i = 0; i < 20; ++i)
        {
            inputs.emplace_back(std::array<double, 4>{
                input[0] + i * 0.1,
                input[1] + i * 0.2,
                input[2] + i * 0.3,
                input[3] + i * 0.4,
            });
        }

        ParameterLib::FunctionParameter<double> param_seq(
            "test", {"x * y + z * t"}, {});
        ParameterLib::FunctionParameter<double> param_par(
            "test", {"x * y + z * t"}, {});

        std::vector<std::vector<double>> seq_results;
        seq_results.reserve(inputs.size());

        for (auto const& arr : inputs)
        {
            double const t = arr[0];
            double const x = arr[1];
            double const y = arr[2];
            double const z = arr[3];
            ParameterLib::SpatialPosition pos;
            pos.setCoordinates(MathLib::Point3d{{x, y, z}});
            seq_results.push_back(param_seq(t, pos));
        }

        std::vector<std::vector<double>> par_results(inputs.size());

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(inputs.size()); ++i)
        {
            auto const& arr = inputs[static_cast<std::size_t>(i)];
            double const t = arr[0];
            double const x = arr[1];
            double const y = arr[2];
            double const z = arr[3];
            ParameterLib::SpatialPosition pos;
            pos.setCoordinates(MathLib::Point3d{{x, y, z}});
            par_results[static_cast<std::size_t>(i)] = param_par(t, pos);
        }

        return seq_results == par_results;
    };

    ac::IntervalGenerator<double> gen(-5.0, 5.0);

    ac::check<std::array<double, 4>>(
        property,
        100,
        ac::make_arbitrary(
            ac::randomTupleGenerator<double, 4, decltype(gen)>{gen}),
        gtest_reporter);
}

TEST_P(FunctionParameterOpenMPTest, IndependentThreadIsolationWithCurves)
{
    auto property = [&](std::vector<std::array<double, 4>> const& inputs)
    {
        ParameterLib::FunctionParameter<double> param(
            "test", {"linear(t) + x + y + z"}, curves_);

        std::vector<std::vector<double>> results(inputs.size());

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(inputs.size()); ++i)
        {
            auto const& arr = inputs[static_cast<std::size_t>(i)];
            double const t = arr[0];
            double const x = arr[1];
            double const y = arr[2];
            double const z = arr[3];
            ParameterLib::SpatialPosition pos;
            pos.setCoordinates(MathLib::Point3d{{x, y, z}});
            results[static_cast<std::size_t>(i)] = param(t, pos);
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

        for (int i = 0; i < static_cast<int>(inputs.size()); ++i)
        {
            auto const& arr = inputs[static_cast<std::size_t>(i)];
            double const t = arr[0];
            double const x = arr[1];
            double const y = arr[2];
            double const z = arr[3];
            double const expected = linear(t) + x + y + z;
            if (std::abs(results[static_cast<std::size_t>(i)][0] - expected) >
                1e-10)
            {
                return false;
            }
        }

        return true;
    };

    ac::check<std::vector<std::array<double, 4>>>(
        property,
        100,
        ac::make_arbitrary(GeneratorForVectorArray4{}),
        gtest_reporter);
}

#ifdef _OPENMP
INSTANTIATE_TEST_SUITE_P(FunctionParameterOpenMP,
                         FunctionParameterOpenMPTest,
                         ::testing::Values(1, 2, 4, 8, 16, 32, 64),
                         [](testing::TestParamInfo<int> const& info)
                         { return "threads_" + std::to_string(info.param); });
#else
INSTANTIATE_TEST_SUITE_P(FunctionParameterOpenMP,
                         FunctionParameterOpenMPTest,
                         ::testing::Values(1),
                         [](testing::TestParamInfo<int> const& info)
                         { return "threads_" + std::to_string(info.param); });
#endif
