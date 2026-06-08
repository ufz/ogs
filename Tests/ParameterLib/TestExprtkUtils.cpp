// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <limits>
#include <map>
#include <memory>
#include <vector>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ParameterLib/ExprtkUtils.h"
#include "ParameterLib/SpatialPosition.h"
#include "Tests/AutoCheckTools.h"

namespace ac = autocheck;

/// Helper to create a simple identity test curve.
static std::unique_ptr<MathLib::PiecewiseLinearInterpolation>
createLinearTestCurve()
{
    std::vector<double> support_points{0.0, 1.0, 2.0};
    std::vector<double> values{0.0, 1.0, 2.0};
    return std::make_unique<MathLib::PiecewiseLinearInterpolation>(
        std::move(support_points), std::move(values));
}

/// Helper to create an alternating test curve.
static std::unique_ptr<MathLib::PiecewiseLinearInterpolation>
createAlternatingTestCurve()
{
    std::vector<double> support_points{0.0, 1.0, 2.0, 3.0};
    std::vector<double> values{0.0, 1.0, 0.0, 1.0};
    return std::make_unique<MathLib::PiecewiseLinearInterpolation>(
        std::move(support_points), std::move(values));
}

struct ExprtkUtilsTest : public ::testing::Test
{
    void SetUp() override
    {
        curves.clear();
        curves["linear"] = createLinearTestCurve();
        curves["alternating"] = createAlternatingTestCurve();
    }

    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;
    ac::gtest_reporter gtest_reporter;
};

// ============================================================================
// CurveWrapper Tests
// ============================================================================

TEST_F(ExprtkUtilsTest, CurveWrapperInSymbolTableProperty)
{
    using namespace ParameterLib;
    // Create symbol table, register curves, and compile expression once
    auto symbol_table = createBaseSymbolTable(false);
    std::map<std::string, CurveWrapper> wrappers;
    registerCurveWrappers(symbol_table, {"linear"}, curves, wrappers);
    auto expressions = compileExpressions(symbol_table, {"linear(t)"});

    auto property = [&symbol_table, &expressions](double const t)
    {
        symbol_table.get_variable("t")->ref() = t;
        return t == expressions[0].value();
    };

    // Test with values in the valid range [0, 2] of the linear curve
    ac::IntervalGenerator<double> generator(0.0, 2.0);
    autocheck::check<double>(property, 100, ac::make_arbitrary(generator),
                             gtest_reporter);
}

// ============================================================================
// SymbolTableCache Tests
// ============================================================================

TEST_F(ExprtkUtilsTest, SymbolTableCacheSetTimeAndPositionWithoutSpatial)
{
    using namespace ParameterLib;
    auto symbol_table = createBaseSymbolTable(false);
    SymbolTableCache cache(symbol_table, false);

    SpatialPosition pos;
    pos.setCoordinates(MathLib::Point3d{{1.0, 2.0, 3.0}});

    // Should not throw even though spatial variables not required
    cache.setTimeAndPosition(5.0, pos);
}

TEST_F(ExprtkUtilsTest, SymbolTableCacheSetTimeAndPositionWithSpatial)
{
    using namespace ParameterLib;
    auto symbol_table = createBaseSymbolTable(true);
    SymbolTableCache cache(symbol_table, true);

    SpatialPosition pos;
    pos.setCoordinates(MathLib::Point3d{{1.0, 2.0, 3.0}});

    // Should succeed - both spatial variables required and provided
    cache.setTimeAndPosition(5.0, pos);

    // Verify values were set by evaluating expressions that use them
    std::vector<std::string> expr_strings{"t", "x", "y", "z"};
    auto expressions = compileExpressions(symbol_table, expr_strings);

    EXPECT_NEAR(5.0, expressions[0].value(),
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(1.0, expressions[1].value(),
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(2.0, expressions[2].value(),
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(3.0, expressions[3].value(),
                std::numeric_limits<double>::epsilon());
}

TEST_F(ExprtkUtilsTest, SymbolTableCacheSetTimeAndPositionMissingSpatialData)
{
    using namespace ParameterLib;
    auto symbol_table = createBaseSymbolTable(true);
    SymbolTableCache cache(symbol_table, true);

    SpatialPosition pos;
    // Do not set coordinates

    // Should throw because spatial variables required but not provided
    EXPECT_THROW(cache.setTimeAndPosition(5.0, pos), std::runtime_error);
}

TEST_F(ExprtkUtilsTest, SymbolTableCacheReinitialize)
{
    using namespace ParameterLib;
    auto symbol_table = createBaseSymbolTable(true);
    SymbolTableCache cache(symbol_table, true);

    // Move the symbol table (simulating relocation)
    auto symbol_table2 = std::move(symbol_table);

    // Reinitialize the cache with the moved symbol table
    cache.reinitialize(symbol_table2);

    SpatialPosition pos;
    pos.setCoordinates(MathLib::Point3d{{2.0, 3.0, 4.0}});
    cache.setTimeAndPosition(7.0, pos);

    std::vector<std::string> expr_strings{"t", "x", "y", "z"};
    auto expressions = compileExpressions(symbol_table2, expr_strings);

    EXPECT_NEAR(7.0, expressions[0].value(),
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(2.0, expressions[1].value(),
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(3.0, expressions[2].value(),
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(4.0, expressions[3].value(),
                std::numeric_limits<double>::epsilon());

    // Reinitialize with same symbol table and verify it works again
    SpatialPosition pos2;
    pos2.setCoordinates(MathLib::Point3d{{4.0, 5.0, 6.0}});
    cache.reinitialize(symbol_table2);
    EXPECT_NO_THROW(cache.setTimeAndPosition(10.0, pos2));
}

// ============================================================================
// hasSpatialVariables Tests
// ============================================================================

TEST_F(ExprtkUtilsTest, HasSpatialVariables)
{
    using namespace ParameterLib;
    // Empty list returns false
    EXPECT_FALSE(hasSpatialVariables(std::vector<std::string>{}));

    // Single spatial coordinates detected
    EXPECT_TRUE(
        hasSpatialVariables(std::vector<std::string>{"t", "x", "temperature"}));
    EXPECT_TRUE(
        hasSpatialVariables(std::vector<std::string>{"t", "y", "pressure"}));
    EXPECT_TRUE(
        hasSpatialVariables(std::vector<std::string>{"t", "z", "velocity"}));

    // Multiple spatial coordinates
    EXPECT_TRUE(
        hasSpatialVariables(std::vector<std::string>{"x", "y", "z", "t"}));

    // No spatial variables returns false
    EXPECT_FALSE(hasSpatialVariables(
        std::vector<std::string>{"t", "temperature", "pressure"}));
}

// ============================================================================
// createBaseSymbolTable Tests
// ============================================================================

TEST_F(ExprtkUtilsTest, CreateBaseSymbolTable)
{
    // Without spatial variables
    {
        auto symbol_table = ParameterLib::createBaseSymbolTable(false);
        EXPECT_NE(nullptr, symbol_table.get_variable("t"));
        EXPECT_EQ(nullptr, symbol_table.get_variable("x"));
        EXPECT_EQ(nullptr, symbol_table.get_variable("y"));
        EXPECT_EQ(nullptr, symbol_table.get_variable("z"));
    }

    // With spatial variables
    {
        auto symbol_table = ParameterLib::createBaseSymbolTable(true);
        EXPECT_NE(nullptr, symbol_table.get_variable("t"));
        EXPECT_NE(nullptr, symbol_table.get_variable("x"));
        EXPECT_NE(nullptr, symbol_table.get_variable("y"));
        EXPECT_NE(nullptr, symbol_table.get_variable("z"));
    }
}

TEST_F(ExprtkUtilsTest, CreateBaseSymbolTableHasConstants)
{
    auto symbol_table = ParameterLib::createBaseSymbolTable(false);
    std::vector<std::string> expr_strings{"pi"};
    auto expressions =
        ParameterLib::compileExpressions(symbol_table, expr_strings);
    double pi_value = expressions[0].value();
    EXPECT_GT(pi_value, 3.0);
    EXPECT_LT(pi_value, 3.3);
}

// ============================================================================
// registerCurveWrappers Tests
// ============================================================================

TEST_F(ExprtkUtilsTest, RegisterCurveWrappers)
{
    using namespace ParameterLib;
    auto symbol_table = createBaseSymbolTable(false);
    std::map<std::string, CurveWrapper> wrappers;
    registerCurveWrappers(symbol_table, {"linear", "alternating"}, curves,
                          wrappers);

    EXPECT_EQ(2, wrappers.size());
    EXPECT_NE(wrappers.end(), wrappers.find("linear"));
    EXPECT_NE(wrappers.end(), wrappers.find("alternating"));

    auto expressions =
        compileExpressions(symbol_table, {"linear(0.5)", "alternating(1.0)"});
    EXPECT_NEAR(0.5, expressions[0].value(),
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(1.0, expressions[1].value(),
                std::numeric_limits<double>::epsilon());
}

// ============================================================================
// compileExpressions Tests
// ============================================================================

TEST_F(ExprtkUtilsTest, CompileExpressions)
{
    using namespace ParameterLib;
    auto symbol_table = createBaseSymbolTable(false);
    symbol_table.create_variable("param");
    symbol_table.get_variable("param")->ref() = 5.0;

    std::map<std::string, CurveWrapper> wrappers;
    registerCurveWrappers(symbol_table, {"linear"}, curves, wrappers);

    // Three expressions: simple arithmetic, with custom variable, with curve
    // function
    auto expressions = compileExpressions(
        symbol_table,
        {"1.0 + 2.0",            // Simple: 3.0
         "param * 2.0",          // Custom variable: 5 * 2 = 10.0
         "linear(1.5) + 0.5"});  // Curve function: 1.5 + 0.5 = 2.0

    EXPECT_EQ(3, expressions.size());
    EXPECT_NEAR(3.0, expressions[0].value(),
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(10.0, expressions[1].value(),
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(2.0, expressions[2].value(),
                std::numeric_limits<double>::epsilon());
}

TEST_F(ExprtkUtilsTest, CompileExpressionsErrors)
{
    using namespace ParameterLib;
    auto symbol_table = createBaseSymbolTable(false);

    // Invalid syntax
    EXPECT_THROW(compileExpressions(symbol_table, {"1.0 +++ 2.0"}),
                 std::runtime_error);

    // Undefined variable
    EXPECT_THROW(compileExpressions(symbol_table, {"undefined_var + 1.0"}),
                 std::runtime_error);

    // Incomplete expression
    EXPECT_THROW(compileExpressions(symbol_table, {"1.0 +"}),
                 std::runtime_error);

    // Unmatched parentheses
    EXPECT_THROW(compileExpressions(symbol_table, {"(1.0 + 2.0"}),
                 std::runtime_error);

    // Mixed valid/invalid expressions
    EXPECT_THROW(
        compileExpressions(symbol_table, {"1.0 + 2.0", "invalid syntax +++"}),
        std::runtime_error);
}

TEST_F(ExprtkUtilsTest, CompileExpressionsRejectsAssignmentToRegisteredSymbol)
{
    using namespace ParameterLib;

    // Assigning to a built-in/registered symbol must fail: it would corrupt
    // the shared per-thread evaluation state.
    {
        auto symbol_table = createBaseSymbolTable(true);
        EXPECT_THROW(compileExpressions(symbol_table, {"x := 5; x"}),
                     std::runtime_error);
    }
    {
        auto symbol_table = createBaseSymbolTable(false);
        EXPECT_THROW(compileExpressions(symbol_table, {"t := 2.0; t"}),
                     std::runtime_error);
    }
    // Assigning to a registered MPL-style variable must also fail.
    {
        auto symbol_table = createBaseSymbolTable(false);
        symbol_table.create_variable("param");
        EXPECT_THROW(compileExpressions(symbol_table, {"param := 1.0; param"}),
                     std::runtime_error);
    }

    // Assignment to a user-defined local 'var' is allowed.
    {
        auto symbol_table = createBaseSymbolTable(false);
        auto expressions =
            compileExpressions(symbol_table, {"var foo := 5; foo + t"});
        symbol_table.get_variable("t")->ref() = 0.0;
        EXPECT_NEAR(5.0, expressions[0].value(),
                    std::numeric_limits<double>::epsilon());
    }
}

// ============================================================================
// collectVariables Tests
// ============================================================================

TEST_F(ExprtkUtilsTest, CollectVariablesMultipleVariables)
{
    std::vector<std::string> expr_strings{"x + y + z + t"};
    auto variables = ParameterLib::collectVariables(expr_strings);

    EXPECT_EQ(4, variables.size());
    EXPECT_NE(std::find(variables.begin(), variables.end(), "x"),
              variables.end());
    EXPECT_NE(std::find(variables.begin(), variables.end(), "y"),
              variables.end());
    EXPECT_NE(std::find(variables.begin(), variables.end(), "z"),
              variables.end());
    EXPECT_NE(std::find(variables.begin(), variables.end(), "t"),
              variables.end());
}

TEST_F(ExprtkUtilsTest, CollectVariablesDeduplication)
{
    std::vector<std::string> expr_strings{"x + x + x", "x + y"};
    auto variables = ParameterLib::collectVariables(expr_strings);

    // x should appear only once despite being used multiple times
    EXPECT_EQ(2, variables.size());
    auto x_count = std::count(variables.begin(), variables.end(), "x");
    EXPECT_EQ(1, x_count);
}

TEST_F(ExprtkUtilsTest, CollectVariablesWithFunctionNames)
{
    // NOTE: collectVariables collects all symbols from expressions, including
    // both variable names and user-defined function (curve) names. This is
    // intentional because exprtk cannot distinguish between variables and
    // user-defined functions in the expression string alone. The distinction
    // is made later by collectUsedCurveNames which filters by checking against
    // the available curves map. This enables detecting which curves are
    // actually used in expressions.
    std::vector<std::string> expr_strings{"linear(x) + alternating(y)"};
    auto variables = ParameterLib::collectVariables(expr_strings);

    // Both variables (x, y) and curve names (linear, alternating) are collected
    EXPECT_NE(std::find(variables.begin(), variables.end(), "linear"),
              variables.end());
    EXPECT_NE(std::find(variables.begin(), variables.end(), "alternating"),
              variables.end());
}

TEST_F(ExprtkUtilsTest, CollectVariablesExcludesConstants)
{
    // The named constants registered by add_constants() (pi, epsilon, inf)
    // must not be reported as variables: downstream code would otherwise try
    // to (re)create them and clash with the constants.
    std::vector<std::string> expr_strings{"pi + epsilon + inf + t"};
    auto variables = ParameterLib::collectVariables(expr_strings);

    EXPECT_EQ(std::find(variables.begin(), variables.end(), "pi"),
              variables.end());
    EXPECT_EQ(std::find(variables.begin(), variables.end(), "epsilon"),
              variables.end());
    EXPECT_EQ(std::find(variables.begin(), variables.end(), "inf"),
              variables.end());
    // t is a normal variable and must still be reported.
    EXPECT_NE(std::find(variables.begin(), variables.end(), "t"),
              variables.end());
}

TEST_F(ExprtkUtilsTest, CollectVariablesEmpty)
{
    std::vector<std::string> expr_strings{};
    auto variables = ParameterLib::collectVariables(expr_strings);

    EXPECT_EQ(0, variables.size());
}

// ============================================================================
// collectUsedCurveNames Tests
// ============================================================================

TEST_F(ExprtkUtilsTest, CollectUsedCurveNamesFiltering)
{
    using namespace ParameterLib;
    // Empty list
    auto used_curves = collectUsedCurveNames({}, curves);
    EXPECT_EQ(0, used_curves.size());

    // No matching curves (non-curve symbols)
    used_curves = collectUsedCurveNames({"x", "y", "z", "t"}, curves);
    EXPECT_EQ(0, used_curves.size());

    // Partial match
    used_curves = collectUsedCurveNames({"x", "linear", "nonexistent"}, curves);
    EXPECT_EQ(1, used_curves.size());
    EXPECT_NE(std::find(used_curves.begin(), used_curves.end(), "linear"),
              used_curves.end());
}

TEST_F(ExprtkUtilsTest, CollectUsedCurveNamesAllMatch)
{
    std::vector<std::string> symbol_names{"linear", "alternating"};
    auto used_curves =
        ParameterLib::collectUsedCurveNames(symbol_names, curves);

    EXPECT_EQ(2, used_curves.size());
    EXPECT_NE(std::find(used_curves.begin(), used_curves.end(), "linear"),
              used_curves.end());
    EXPECT_NE(std::find(used_curves.begin(), used_curves.end(), "alternating"),
              used_curves.end());
}

// ============================================================================
// compileExpressions Error Handling
// ============================================================================

TEST_F(ExprtkUtilsTest, CompileExpressionsEmptyList)
{
    auto symbol_table = ParameterLib::createBaseSymbolTable(false);
    std::vector<std::string> expr_strings{};

    auto expressions =
        ParameterLib::compileExpressions(symbol_table, expr_strings);
    EXPECT_EQ(0, expressions.size());
}

// ============================================================================
// SymbolTableCache Pointer Safety & State
// ============================================================================

TEST_F(ExprtkUtilsTest, SymbolTableCacheDefaultConstructor)
{
    // Test default constructor - should create inactive cache without throwing
    // Default cache should work for lazy initialization patterns
    EXPECT_NO_THROW(ParameterLib::SymbolTableCache());
}

TEST_F(ExprtkUtilsTest, SymbolTableCacheWithSpatialThenWithoutConsistency)
{
    using namespace ParameterLib;
    // Test that switching between spatial and non-spatial is safe
    auto symbol_table_no_spatial = createBaseSymbolTable(false);
    SymbolTableCache cache_no_spatial(symbol_table_no_spatial, false);

    SpatialPosition pos;
    pos.setCoordinates(MathLib::Point3d{{1.0, 2.0, 3.0}});

    // Should not throw - spatial data is ignored when not required
    EXPECT_NO_THROW(cache_no_spatial.setTimeAndPosition(5.0, pos));
}

// ============================================================================
// registerCurveWrappers Error Handling
// ============================================================================

TEST_F(ExprtkUtilsTest, RegisterCurveWrappersErrors)
{
    using namespace ParameterLib;
    auto symbol_table = createBaseSymbolTable(false);
    std::map<std::string, CurveWrapper> wrappers;

    // Nonexistent curve
    EXPECT_THROW(
        registerCurveWrappers(symbol_table, {"nonexistent"}, curves, wrappers),
        std::out_of_range);

    // Mix of valid and invalid curve names
    EXPECT_THROW(registerCurveWrappers(symbol_table, {"linear", "nonexistent"},
                                       curves, wrappers),
                 std::out_of_range);
}

TEST_F(ExprtkUtilsTest, RegisterCurveWrappersRejectsReservedName)
{
    using namespace ParameterLib;
    auto symbol_table = createBaseSymbolTable(false);
    std::map<std::string, CurveWrapper> wrappers;

    // A curve named like an exprtk built-in (e.g. "sin") cannot be registered:
    // exprtk would silently keep the built-in and the curve would be ignored.
    // Such a name must be rejected outright.
    curves["sin"] = createLinearTestCurve();
    EXPECT_THROW(registerCurveWrappers(symbol_table, {"sin"}, curves, wrappers),
                 std::runtime_error);
}

// ============================================================================
// State Management & Cache Reuse
// ============================================================================

TEST_F(ExprtkUtilsTest, SymbolTableCacheCrossTimeStepReuse)
{
    using namespace ParameterLib;
    // Verify cache works correctly across multiple time steps
    auto symbol_table = createBaseSymbolTable(true);
    SymbolTableCache cache(symbol_table, true);

    // Compile expression once
    std::vector<std::string> expr_strings{"x + y + z + t"};
    auto expressions = compileExpressions(symbol_table, expr_strings);

    // Evaluate at multiple time steps
    for (int step = 0; step < 5; ++step)
    {
        SpatialPosition pos;
        pos.setCoordinates(MathLib::Point3d{{static_cast<double>(step),
                                             static_cast<double>(step + 1),
                                             static_cast<double>(step + 2)}});
        cache.setTimeAndPosition(static_cast<double>(step * 10), pos);

        double expected = step + (step + 1) + (step + 2) + (step * 10);
        EXPECT_NEAR(expected, expressions[0].value(),
                    std::numeric_limits<double>::epsilon());
    }
}

TEST_F(ExprtkUtilsTest, CollectVariablesComplexExpressions)
{
    // Test with nested function calls and complex expressions
    // Note: collectVariables only collects user-defined variables and curves,
    // not built-in functions (sin, cos, sqrt are language primitives)
    std::vector<std::string> expr_strings{
        "sin(x) * cos(y) + linear(z)", "sqrt(t) + alternating(x + y)",
        "(linear(z) - alternating(t)) / (x * y)"};

    auto variables = ParameterLib::collectVariables(expr_strings);

    // Should collect spatial/temporal variables
    EXPECT_NE(std::find(variables.begin(), variables.end(), "x"),
              variables.end());
    EXPECT_NE(std::find(variables.begin(), variables.end(), "y"),
              variables.end());
    EXPECT_NE(std::find(variables.begin(), variables.end(), "z"),
              variables.end());
    EXPECT_NE(std::find(variables.begin(), variables.end(), "t"),
              variables.end());

    // Should collect user-defined curve names
    EXPECT_NE(std::find(variables.begin(), variables.end(), "linear"),
              variables.end());
    EXPECT_NE(std::find(variables.begin(), variables.end(), "alternating"),
              variables.end());

    // Built-in functions (sin, cos, sqrt) are not collected by exprtk
    // as they are language primitives, not variables
}

TEST_F(ExprtkUtilsTest, CollectUsedCurveNamesComplexExpressions)
{
    // Test with the complex expressions from above
    std::vector<std::string> expr_strings{
        "sin(x) * cos(y) + linear(z)", "sqrt(t) + alternating(x + y)",
        "(linear(z) - alternating(t)) / (x * y)"};

    auto variables = ParameterLib::collectVariables(expr_strings);
    auto curve_names = ParameterLib::collectUsedCurveNames(variables, curves);

    // Should identify only the curves, not built-in functions
    EXPECT_EQ(2, curve_names.size());
    EXPECT_NE(std::find(curve_names.begin(), curve_names.end(), "linear"),
              curve_names.end());
    EXPECT_NE(std::find(curve_names.begin(), curve_names.end(), "alternating"),
              curve_names.end());

    // Should NOT include built-in functions
    EXPECT_EQ(std::find(curve_names.begin(), curve_names.end(), "sin"),
              curve_names.end());
    EXPECT_EQ(std::find(curve_names.begin(), curve_names.end(), "cos"),
              curve_names.end());
}

TEST_F(ExprtkUtilsTest, SymbolTableCacheErrorsAreDetected)
{
    using namespace ParameterLib;
    // Verify that cache detects and reports errors when spatial coords are
    // required but missing. Note: implementation may update time before
    // throwing error on missing coordinates (non-atomic operation).
    auto symbol_table = createBaseSymbolTable(true);
    SymbolTableCache cache(symbol_table, true);

    SpatialPosition valid_pos;
    valid_pos.setCoordinates(MathLib::Point3d{{1.0, 2.0, 3.0}});
    cache.setTimeAndPosition(5.0, valid_pos);

    // Verify values were set
    std::vector<std::string> expr_strings{"t", "x", "y", "z"};
    auto expressions = compileExpressions(symbol_table, expr_strings);
    EXPECT_NEAR(5.0, expressions[0].value(),
                std::numeric_limits<double>::epsilon());

    // Try to set without coordinates (should throw)
    SpatialPosition invalid_pos;
    EXPECT_THROW(cache.setTimeAndPosition(10.0, invalid_pos),
                 std::runtime_error);

    // Spatial coordinates should remain at their previous values
    expressions = compileExpressions(symbol_table, expr_strings);
    EXPECT_NEAR(1.0, expressions[1].value(),
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(2.0, expressions[2].value(),
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(3.0, expressions[3].value(),
                std::numeric_limits<double>::epsilon());
}

// ============================================================================
// Integration Tests
// ============================================================================

TEST_F(ExprtkUtilsTest, FullWorkflowWithoutSpatialVariables)
{
    using namespace ParameterLib;
    // Create symbol table
    auto symbol_table = createBaseSymbolTable(false);

    // Add a custom variable. In OGS this stands in for an MPL variable: in
    // MaterialLib/MPL/Properties/Function.cpp such symbols are bound to
    // VariableArray fields via symbol_table.add_variable / add_vector. Here we
    // create and set it directly to simulate that binding.
    symbol_table.create_variable("custom");
    symbol_table.get_variable("custom")->ref() = 3.0;

    // Register curves
    std::map<std::string, CurveWrapper> wrappers;
    registerCurveWrappers(symbol_table, {"linear"}, curves, wrappers);

    // Compile expressions
    std::vector<std::string> expr_strings{"custom * linear(2.0) + t"};
    auto expressions = compileExpressions(symbol_table, expr_strings);

    // Set time and evaluate
    symbol_table.get_variable("t")->ref() = 1.0;
    double result = expressions[0].value();

    // custom (3) * linear(2) (2) + t (1) = 6 + 1 = 7
    EXPECT_NEAR(7.0, result, std::numeric_limits<double>::epsilon());
}

TEST_F(ExprtkUtilsTest, FullWorkflowWithSpatialVariables)
{
    using namespace ParameterLib;
    // Create symbol table with spatial variables
    auto symbol_table = createBaseSymbolTable(true);

    // Register curves
    std::map<std::string, CurveWrapper> wrappers;
    registerCurveWrappers(symbol_table, {"alternating"}, curves, wrappers);

    // Create cache for efficient evaluation
    SymbolTableCache cache(symbol_table, true);

    // Compile expressions
    std::vector<std::string> expr_strings{"x + y + z + alternating(t)"};
    auto expressions = compileExpressions(symbol_table, expr_strings);

    // Set position and time
    SpatialPosition pos;
    pos.setCoordinates(MathLib::Point3d{{1.0, 2.0, 3.0}});
    cache.setTimeAndPosition(1.0, pos);

    double result = expressions[0].value();

    // x (1) + y (2) + z (3) + alternating(1) (1) = 7
    EXPECT_NEAR(7.0, result, std::numeric_limits<double>::epsilon());
}

TEST_F(ExprtkUtilsTest, CollectVariablesAndCurvesIntegration)
{
    // Simulate real-world expression analysis
    std::vector<std::string> expr_strings{"linear(t) * x + y",
                                          "alternating(z) + 1.0"};

    // Collect all variables
    auto variables = ParameterLib::collectVariables(expr_strings);

    // Determine which variables are curves
    auto curve_names = ParameterLib::collectUsedCurveNames(variables, curves);

    // Verify correct identification
    EXPECT_NE(std::find(variables.begin(), variables.end(), "linear"),
              variables.end());
    EXPECT_NE(std::find(variables.begin(), variables.end(), "alternating"),
              variables.end());
    EXPECT_NE(std::find(variables.begin(), variables.end(), "x"),
              variables.end());
    EXPECT_NE(std::find(variables.begin(), variables.end(), "y"),
              variables.end());
    EXPECT_NE(std::find(variables.begin(), variables.end(), "z"),
              variables.end());

    EXPECT_EQ(2, curve_names.size());
    EXPECT_NE(std::find(curve_names.begin(), curve_names.end(), "linear"),
              curve_names.end());
    EXPECT_NE(std::find(curve_names.begin(), curve_names.end(), "alternating"),
              curve_names.end());
}
