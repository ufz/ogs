// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "ExprtkUtils.h"
#include "SpatialPosition.h"

namespace ParameterLib
{

/// Shared implementation for mathematical function evaluation using exprtk.
/// Provides per-thread storage for thread-safe evaluation under OpenMP.
class FunctionEvaluation
{
public:
    using Expression = exprtk::expression<double>;
    using SymbolTable = exprtk::symbol_table<double>;

    /// \param num_threads  Number of threads for per-thread storage.
    /// \param variables    Variable names to register in the symbol table.
    ///                     Spatial variables (x, y, z) are created only when
    ///                     they appear here.
    /// \param expression_strings  Expression strings to compile.
    /// \param curves       Named piecewise-linear curves available to
    ///                     expressions.
    FunctionEvaluation(
        int num_threads,
        std::vector<std::string> const& variables,
        std::vector<std::string> const& expression_strings,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves);

    /// Evaluate all value expressions at the given spatial position and time.
    ///
    /// Although marked \c const, this method writes to per-thread symbol-table
    /// slots through cached pointers. Each OpenMP thread writes exclusively to
    /// its own slot, so concurrent calls from different threads do not race.
    /// The mathematical function itself is immutable after construction.
    std::vector<double> evaluate(SpatialPosition const& pos, double t) const;

    /// Evaluate all value expressions into a caller-provided output buffer.
    ///
    /// No heap allocation is performed. Thread-safety is identical to the
    /// vector-returning overload.
    /// \param pos    Spatial position for evaluation.
    /// \param t      Time for evaluation.
    /// \param result Must have size equal to getNumberOfComponents().
    void evaluate(SpatialPosition const& pos, double t,
                  std::span<double> result) const;

    /// Number of components (number of compiled value expressions).
    int getNumberOfComponents() const;

    bool isTimeDependent() const { return true; }

private:
    /// Symbol table storage and compiled expressions for one OpenMP thread.
    /// Each thread receives its own instance so that concurrent evaluations
    /// do not race on the exprtk symbol table variables.
    struct PerThreadData
    {
        PerThreadData(
            std::vector<std::string> const& variables,
            bool spatial_position_is_required,
            std::vector<std::string> const& used_curve_names,
            std::map<
                std::string,
                std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
                curves,
            std::vector<std::string> const& expression_strings);

        PerThreadData(PerThreadData&& other);
        PerThreadData& operator=(PerThreadData&& other) noexcept = delete;
        PerThreadData(PerThreadData const&) = delete;
        PerThreadData& operator=(PerThreadData const&) = delete;

        /// Curve wrappers owned by this thread; must outlive symbol_table.
        std::map<std::string, CurveWrapper> curve_wrappers;
        /// Symbol table holding variable storage. Must be destroyed after
        /// value_expressions (exprtk requirement) and after curve_wrappers.
        SymbolTable symbol_table;
        /// Compiled value expressions. Must be destroyed before symbol_table.
        std::vector<Expression> value_expressions;
        /// Cached pointers to t, x, y, z inside symbol_table.
        SymbolTableCache symbol_table_cache;
    };

    /// Per-thread evaluation context; indexed by omp_get_thread_num().
    std::vector<PerThreadData> per_thread_data_;
};

}  // namespace ParameterLib
