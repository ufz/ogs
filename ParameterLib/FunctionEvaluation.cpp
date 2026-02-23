// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "FunctionEvaluation.h"

#include <cassert>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "BaseLib/Error.h"
#include "ExprtkUtils.h"

namespace ParameterLib
{

FunctionEvaluation::FunctionEvaluation(
    int num_threads,
    std::vector<std::string> const& variables,
    std::vector<std::string> const& expression_strings,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    auto const spatial_position_is_required = hasSpatialVariables(variables);

    auto const expression_symbol_names = collectVariables(expression_strings);
    auto const used_curve_names =
        collectUsedCurveNames(expression_symbol_names, curves);

    per_thread_data_.reserve(num_threads);
    for (int thread_id = 0; thread_id < num_threads; ++thread_id)
    {
        per_thread_data_.emplace_back(variables, spatial_position_is_required,
                                      used_curve_names, curves,
                                      expression_strings);
    }
}

static FunctionEvaluation::SymbolTable createSymbolTable(
    std::vector<std::string> const& variables,
    bool const spatial_position_is_required,
    std::vector<std::string> const& used_curve_names,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::map<std::string, CurveWrapper>& curve_wrappers)
{
    auto symbol_table = createBaseSymbolTable(spatial_position_is_required);

    // Create additional (non-standard) variables.
    for (auto const& v : variables)
    {
        if (!isBuiltinSymbol(v))
        {
            symbol_table.create_variable(v);
        }
    }

    registerCurveWrappers(symbol_table, used_curve_names, curves,
                          curve_wrappers);

    return symbol_table;
}

FunctionEvaluation::PerThreadData::PerThreadData(
    std::vector<std::string> const& variables,
    bool const spatial_position_is_required,
    std::vector<std::string> const& used_curve_names,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::vector<std::string> const& expression_strings)
    : symbol_table(createSymbolTable(variables, spatial_position_is_required,
                                     used_curve_names, curves, curve_wrappers)),
      value_expressions(compileExpressions(symbol_table, expression_strings)),
      symbol_table_cache(symbol_table, spatial_position_is_required)
{
}

FunctionEvaluation::PerThreadData::PerThreadData(PerThreadData&& /*other*/)
{
    OGS_FATAL(
        "Function parameter: The internal PerThreadData is not "
        "move-constructible.");
}

void FunctionEvaluation::evaluate(SpatialPosition const& pos, double const t,
                                  std::span<double> result) const
{
#ifdef _OPENMP
    int const thread_id = omp_get_thread_num();
#else
    int const thread_id = 0;
#endif

    if (thread_id >= static_cast<int>(per_thread_data_.size()))
    {
        OGS_FATAL("Thread id {:d} exceeds allocated threads {:d}.", thread_id,
                  per_thread_data_.size());
    }

    auto const& thread_data = per_thread_data_[thread_id];
    auto const& expressions = thread_data.value_expressions;
    assert(result.size() == expressions.size());

    auto const& cache = thread_data.symbol_table_cache;

    // Update time and spatial position.
    cache.setTimeAndPosition(t, pos);

    // Evaluate expressions.
    for (unsigned i = 0; i < expressions.size(); ++i)
    {
        result[i] = expressions[i].value();
    }
}

std::vector<double> FunctionEvaluation::evaluate(SpatialPosition const& pos,
                                                 double const t) const
{
    std::vector<double> result(
        per_thread_data_.front().value_expressions.size());
    evaluate(pos, t, result);
    return result;
}

int FunctionEvaluation::getNumberOfComponents() const
{
    return static_cast<int>(per_thread_data_.front().value_expressions.size());
}

}  // namespace ParameterLib
