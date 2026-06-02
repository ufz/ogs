// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ExprtkUtils.h"

#include <cassert>

#include "BaseLib/Algorithm.h"

namespace ParameterLib
{

void SymbolTableCache::setTimeAndPosition(double const t,
                                          SpatialPosition const& pos) const
{
    // t_ptr is null only on a default-constructed cache, which exists solely
    // to satisfy compilation of PerThreadData's (forbidden, OGS_FATAL) move
    // ctor and is never evaluated against.
    assert(t_ptr != nullptr);
    *t_ptr = t;
    // x_ptr is non-null iff the constructor was called with
    // spatial_position_is_required=true; use it as the sentinel.
    if (x_ptr != nullptr)
    {
        // Invariant: all three spatial pointers are set together.
        assert(y_ptr != nullptr && z_ptr != nullptr);
        if (!pos.getCoordinates())
        {
            OGS_FATAL(
                "An expression requires spatial position (x, y, z), "
                "but no coordinates have been set.");
        }
        auto const coords = pos.getCoordinates().value();
        *x_ptr = coords[0];
        *y_ptr = coords[1];
        *z_ptr = coords[2];
    }
}

bool isBuiltinSymbol(std::string_view const name)
{
    return name == "t" || name == "x" || name == "y" || name == "z";
}

bool hasSpatialVariables(std::vector<std::string> const& variables)
{
    return std::ranges::any_of(variables, [](std::string const& v)
                               { return v == "x" || v == "y" || v == "z"; });
}

exprtk::symbol_table<double> createBaseSymbolTable(
    bool spatial_position_is_required)
{
    exprtk::symbol_table<double> symbol_table;
    symbol_table.add_constants();
    symbol_table.create_variable("t");
    if (spatial_position_is_required)
    {
        symbol_table.create_variable("x");
        symbol_table.create_variable("y");
        symbol_table.create_variable("z");
    }
    return symbol_table;
}

void registerCurveWrappers(
    exprtk::symbol_table<double>& symbol_table,
    std::vector<std::string> const& curve_names,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::map<std::string, CurveWrapper>& curve_wrappers)
{
    for (auto const& name : curve_names)
    {
        curve_wrappers.emplace(name, CurveWrapper(*curves.at(name)));
    }
    for (auto& [name, wrapper] : curve_wrappers)
    {
        symbol_table.add_function(name, wrapper);
    }
}

SymbolTableCache::SymbolTableCache(exprtk::symbol_table<double>& symbol_table,
                                   bool const spatial_position_is_required)
    : spatial_position_is_required_(spatial_position_is_required)
{
    reinitialize(symbol_table);
}

void SymbolTableCache::reinitialize(exprtk::symbol_table<double>& symbol_table)
{
    t_ptr = &(symbol_table.get_variable("t")->ref());
    if (spatial_position_is_required_)
    {
        x_ptr = &(symbol_table.get_variable("x")->ref());
        y_ptr = &(symbol_table.get_variable("y")->ref());
        z_ptr = &(symbol_table.get_variable("z")->ref());
    }
    else
    {
        x_ptr = nullptr;
        y_ptr = nullptr;
        z_ptr = nullptr;
    }
}

template std::vector<exprtk::expression<double>> compileExpressions<double>(
    exprtk::symbol_table<double>&, std::vector<std::string> const&);

std::vector<std::string> collectVariables(
    std::vector<std::string> const& expression_strings)
{
    std::vector<std::string> expression_symbol_names;
    for (auto const& expr : expression_strings)
    {
        if (!exprtk::collect_variables(expr, expression_symbol_names))
        {
            OGS_FATAL("Collecting variables from expression '{}' didn't work.",
                      expr);
        }
    }
    BaseLib::makeVectorUnique(expression_symbol_names);

    // exprtk::collect_variables reports the named constants added by
    // add_constants() (pi, epsilon, inf) as variables. Drop them so that
    // downstream code does not try to (re)create them as variables and clash
    // with the constants. A base symbol table is used so the set of constants
    // is queried from exprtk itself rather than hard-coded.
    auto const constants = createBaseSymbolTable(false);
    std::erase_if(expression_symbol_names,
                  [&constants](std::string const& name)
                  { return constants.is_constant_node(name); });

    return expression_symbol_names;
}

std::vector<std::string> collectUsedCurveNames(
    std::vector<std::string> const& expression_symbol_names,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    std::vector<std::string> used;
    for (auto const& v : expression_symbol_names)
    {
        if (curves.contains(v))
        {
            used.push_back(v);
        }
    }
    return used;
}

}  // namespace ParameterLib
