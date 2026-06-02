// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <exprtk.hpp>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "BaseLib/Error.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "SpatialPosition.h"

namespace ParameterLib
{

/// Wrapper for curve interpolation to be used with exprtk.
/// Converts MathLib::PiecewiseLinearInterpolation to an exprtk function
/// interface.
struct CurveWrapper : public exprtk::ifunction<double>
{
    explicit CurveWrapper(MathLib::PiecewiseLinearInterpolation const& curve)
        : exprtk::ifunction<double>(1), curve_(curve)
    {
        exprtk::disable_has_side_effects(*this);
    }

    double operator()(double const& t) override { return curve_.getValue(t); }

private:
    MathLib::PiecewiseLinearInterpolation const& curve_;
};

/// Cached pointers to the t, x, y, z symbol-table variables per thread.
/// Eliminates string-based lookups during evaluation.
struct SymbolTableCache
{
    SymbolTableCache() = default;

    /// Caches pointers to t and, if \p spatial_position_is_required, to x, y,
    /// z from \p symbol_table. The invariant x_ptr/y_ptr/z_ptr are non-null
    /// iff spatial_position_is_required is true is established here.
    SymbolTableCache(exprtk::symbol_table<double>& symbol_table,
                     bool spatial_position_is_required);

    /// Re-caches pointers from a (possibly moved) symbol_table.
    void reinitialize(exprtk::symbol_table<double>& symbol_table);

    /// Sets t and, when spatial variables were registered
    /// (x/y/z_ptr != nullptr), also sets x, y, z from \p pos.
    void setTimeAndPosition(double t, SpatialPosition const& pos) const;

private:
    bool spatial_position_is_required_ = false;
    double* t_ptr = nullptr;
    double* x_ptr = nullptr;
    double* y_ptr = nullptr;
    double* z_ptr = nullptr;
};

/// Compiles expression strings into exprtk expressions using the given symbol
/// table.
/// Assignment (`:=`, `+=`, ...) to user-defined local `var`s is allowed, but
/// assignment to a symbol that is already registered in \p symbol_table
/// (the built-in t/x/y/z, MPL variables, curves or constants) is rejected:
/// such a write would corrupt the shared per-thread evaluation scratch. This
/// is enforced by enabling exprtk's dependent-entity collector and checking,
/// after each compile, whether any assignment target resolves to a registered
/// symbol (local `var` assignments are not reported by the collector).
template <typename T>
std::vector<exprtk::expression<T>> compileExpressions(
    exprtk::symbol_table<T>& symbol_table,
    std::vector<std::string> const& string_expressions)
{
    using settings_t = typename exprtk::parser<T>::settings_store;
    exprtk::parser<T> parser(settings_t::default_compile_all_opts +
                             settings_t::e_collect_assings);

    std::vector<exprtk::expression<T>> expressions(string_expressions.size());
    for (unsigned i = 0; i < string_expressions.size(); ++i)
    {
        expressions[i].register_symbol_table(symbol_table);
        if (!parser.compile(string_expressions[i], expressions[i]))
        {
            OGS_FATAL("Error: {:s}\tExpression: {:s}\n", parser.error(),
                      string_expressions[i]);
        }

        std::vector<
            typename exprtk::parser<T>::dependent_entity_collector::symbol_t>
            assignments;
        parser.dec().assignment_symbols(assignments);
        if (!assignments.empty())
        {
            std::string names;
            for (auto const& [name, type] : assignments)
            {
                if (!names.empty())
                {
                    names += ", ";
                }
                names += name;
            }
            OGS_FATAL(
                "Expression '{:s}' assigns to the already-defined symbol(s) "
                "'{:s}'. Assignment to the built-in variables (t, x, y, z), to "
                "MPL variables, curves or constants is not allowed because it "
                "would corrupt the shared evaluation state. Assignment to "
                "user-defined local 'var's is permitted.",
                string_expressions[i], names);
        }
    }
    return expressions;
}

/// Returns true if \p name is one of the built-in symbol names registered by
/// createBaseSymbolTable ("t", "x", "y", "z").
bool isBuiltinSymbol(std::string_view name);

/// Returns true if \p variables contains any of the spatial variable names
/// "x", "y", or "z".
bool hasSpatialVariables(std::vector<std::string> const& variables);

/// Creates a symbol table with time variable ("t") and optionally spatial
/// variables ("x", "y", "z"). Does not add any other variables or curves.
exprtk::symbol_table<double> createBaseSymbolTable(
    bool spatial_position_is_required);

/// Registers CurveWrapper objects in the symbol table for the given curve
/// names.
/// \param symbol_table  The symbol table to register curves into.
/// \param curve_names   Names of curves to register.
/// \param curves        Map of all available curves.
/// \param curve_wrappers Is filled with CurveWrapper objects and must
///                       outlive the returned symbol table.
void registerCurveWrappers(
    exprtk::symbol_table<double>& symbol_table,
    std::vector<std::string> const& curve_names,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::map<std::string, CurveWrapper>& curve_wrappers);

extern template std::vector<exprtk::expression<double>>
compileExpressions<double>(exprtk::symbol_table<double>&,
                           std::vector<std::string> const&);

/// Collects and deduplicates all variable/function identifiers appearing in
/// the given expression strings by calling exprtk::collect_variables on each.
std::vector<std::string> collectVariables(
    std::vector<std::string> const& expression_strings);

/// Returns the subset of \p expression_symbol_names that are keys in curves.
std::vector<std::string> collectUsedCurveNames(
    std::vector<std::string> const& expression_symbol_names,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);

}  // namespace ParameterLib
