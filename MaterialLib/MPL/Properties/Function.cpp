/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MaterialLib/MPL/Properties/Function.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <exprtk.hpp>
#include <numeric>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>
#include <ranges>
#include <typeinfo>
#include <unordered_set>

#include "BaseLib/Algorithm.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/VectorizedTensor.h"

namespace MaterialPropertyLib
{
struct CurveWrapper : public exprtk::ifunction<double>
{
    explicit CurveWrapper(const MathLib::PiecewiseLinearInterpolation& curve)
        : exprtk::ifunction<double>(1), _curve(curve)
    {
        exprtk::disable_has_side_effects(*this);
    }

    double operator()(const double& t) override { return _curve.getValue(t); }

private:
    const MathLib::PiecewiseLinearInterpolation& _curve;
};
// Passing symbol table by reference as required by register_symbol_table()
// call.
template <typename T>
static std::vector<exprtk::expression<T>> compileExpressions(
    exprtk::symbol_table<T>& symbol_table,
    std::vector<std::string> const& string_expressions)
{
    exprtk::parser<T> parser;

    std::vector<exprtk::expression<T>> expressions(string_expressions.size());
    for (unsigned i = 0; i < string_expressions.size(); ++i)
    {
        expressions[i].register_symbol_table(symbol_table);
        if (!parser.compile(string_expressions[i], expressions[i]))
        {
            OGS_FATAL("Error: {:s}\tExpression: {:s}\n",
                      parser.error(),
                      string_expressions[i]);
        }
    }
    return expressions;
}

template <int D>
class Function::Implementation
{
public:
    using Expression = exprtk::expression<double>;

public:
    Implementation(
        int num_threads,
        std::vector<std::string> const& variables,
        std::vector<std::string> const& value_string_expressions,
        std::vector<std::pair<std::string, std::vector<std::string>>> const&
            dvalue_string_expressions,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves);

private:
    /// Create symbol table for given variables and populates the variable_array
    /// as needed.
    exprtk::symbol_table<double> createSymbolTable(
        std::vector<std::string> const& variables,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves,
        VariableArray& variable_array);

    /// Create symbol tables for all threads.
    std::vector<exprtk::symbol_table<double>> createSymbolTables(
        int num_threads,
        std::vector<std::string> const& variables,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves);

public:
    /// Value expressions per thread.
    /// Multiple expressions are representing vector-valued functions.
    std::vector<std::vector<Expression>> value_expressions;

    /// Derivative expressions with respect to the variable per thread.
    /// Multiple expressions are representing vector-valued functions.
    std::vector<std::vector<std::pair<Variable, std::vector<Expression>>>>
        dvalue_expressions;

    /// Stores values for evaluation of vectorial quantities per thread. Needed
    /// for constant pointers for exprtk.
    std::vector<VariableArray> variable_arrays;

    std::map<std::string, CurveWrapper> _curve_wrappers;

    bool spatial_position_is_required = false;
};

template <int D>
exprtk::symbol_table<double> Function::Implementation<D>::createSymbolTable(
    std::vector<std::string> const& variables,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    VariableArray& variable_array)
{
    exprtk::symbol_table<double> symbol_table;
    symbol_table.add_constants();

    std::unordered_set<std::string> curve_names;
    for (auto const& curve : curves)
    {
        curve_names.insert(curve.first);
    }

    std::unordered_set<std::string> used_curves;

    symbol_table.create_variable("t");

    for (auto const& v : variables)
    {
        if (v == "t")
        {
            continue;
        }
        else if (v == "x" || v == "y" || v == "z")
        {
            symbol_table.create_variable("x");
            symbol_table.create_variable("y");
            symbol_table.create_variable("z");
            spatial_position_is_required = true;
        }
        else if (curve_names.contains(v))
        {
            used_curves.insert(v);
        }
        else
        {
            auto add_scalar = [&v, &symbol_table](double& value)
            { symbol_table.add_variable(v, value); };

            auto add_vector =
                [&v, &symbol_table](double* ptr, std::size_t const size)
            { symbol_table.add_vector(v, ptr, size); };

            auto add_any_variable = BaseLib::Overloaded{
                [&add_scalar](VariableArray::Scalar* address)
                { add_scalar(*address); },
                [&add_vector](VariableArray::KelvinVector* address)
                {
                    auto constexpr size =
                        MathLib::KelvinVector::kelvin_vector_dimensions(D);
                    auto& result = address->template emplace<
                        Eigen::Matrix<double, size, 1>>();
                    add_vector(result.data(), size);
                },
                [&add_vector](VariableArray::DeformationGradient* address)
                {
                    auto constexpr size = MathLib::VectorizedTensor::size(D);
                    auto& result = address->template emplace<
                        Eigen::Matrix<double, size, 1>>();
                    add_vector(result.data(), size);
                }};

            Variable const variable = convertStringToVariable(v);
            variable_array.visitVariable(add_any_variable, variable);
        }
    }

    for (const auto& name : used_curves)
    {
        const auto& curve_ptr = curves.at(name);
        _curve_wrappers.emplace(name, CurveWrapper(*curve_ptr));
    }
    for (auto& [name, wrapper] : _curve_wrappers)
    {
        symbol_table.add_function(name, wrapper);
    }
    return symbol_table;
}
template <int D>
std::vector<exprtk::symbol_table<double>>
Function::Implementation<D>::createSymbolTables(
    int num_threads,
    std::vector<std::string> const& variables,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    std::vector<exprtk::symbol_table<double>> symbol_tables(num_threads);

    for (int thread_id = 0; thread_id < num_threads; ++thread_id)
    {
        symbol_tables[thread_id] =
            createSymbolTable(variables, curves, variable_arrays[thread_id]);
    }

    return symbol_tables;
}

template <int D>
Function::Implementation<D>::Implementation(
    int num_threads,
    std::vector<std::string> const& variables,
    std::vector<std::string> const& value_string_expressions,
    std::vector<std::pair<std::string, std::vector<std::string>>> const&
        dvalue_string_expressions,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    variable_arrays.resize(num_threads);
    value_expressions.resize(num_threads);
    dvalue_expressions.resize(num_threads);

    auto symbol_tables = createSymbolTables(num_threads, variables, curves);

    // value expressions per thread.
    for (int thread_id = 0; thread_id < num_threads; ++thread_id)
    {
        value_expressions[thread_id] = compileExpressions(
            symbol_tables[thread_id], value_string_expressions);
    }

    // dValue expressions per thread.
    for (int thread_id = 0; thread_id < num_threads; ++thread_id)
    {
        for (auto const& [variable_name, string_expressions] :
             dvalue_string_expressions)
        {
            dvalue_expressions[thread_id].emplace_back(
                convertStringToVariable(variable_name),
                compileExpressions(symbol_tables[thread_id],
                                   string_expressions));
        }
    }
}

static void updateVariableArrayValues(std::vector<Variable> const& variables,
                                      VariableArray const& new_variable_array,
                                      VariableArray& variable_array)
{
    for (auto const& variable : variables)
    {
        auto assign_scalar =
            [&variable, &new_variable_array](VariableArray::Scalar* address)
        {
            double const value = *std::get<VariableArray::Scalar const*>(
                new_variable_array.address_of(variable));

            if (std::isnan(value))
            {
                OGS_FATAL(
                    "Function property: Scalar variable '{:s}' is not "
                    "initialized.",
                    variable_enum_to_string[static_cast<int>(variable)]);
            }

            *address = value;
        };
        auto assign_kelvin_vector = [&variable, &new_variable_array](
                                        VariableArray::KelvinVector* address)
        {
            auto assign_value = [&destination = *address,
                                 &variable]<typename S>(S const& source)
            {
                if constexpr (std::is_same_v<S, std::monostate>)
                {
                    OGS_FATAL(
                        "Function property: Kelvin vector variable '{:s}' is "
                        "not initialized.",
                        variable_enum_to_string[static_cast<int>(variable)]);
                }
                else
                {
                    if (std::holds_alternative<S>(destination))
                    {
                        std::get<S>(destination) = MathLib::KelvinVector::
                            kelvinVectorToSymmetricTensor(source);
                    }
                    else
                    {
                        OGS_FATAL(
                            "Function property: Mismatch of Kelvin vector "
                            "sizes for variable {:s}.",
                            variable_enum_to_string[static_cast<int>(
                                variable)]);
                    }
                }
            };
            std::visit(assign_value,
                       *std::get<VariableArray::KelvinVector const*>(
                           new_variable_array.address_of(variable)));
        };
        auto assign_deformation_gradient =
            [&variable,
             &new_variable_array](VariableArray::DeformationGradient* address)
        {
            auto assign_value = [&destination = *address,
                                 &variable]<typename S>(S const& source)
            {
                if constexpr (std::is_same_v<S, std::monostate>)
                {
                    OGS_FATAL(
                        "Function property: Vectorized tensor variable '{:s}' "
                        "is not initialized.",
                        variable_enum_to_string[static_cast<int>(variable)]);
                }
                else
                {
                    if (std::holds_alternative<S>(destination))
                    {
                        std::get<S>(destination) = source;
                    }
                    else
                    {
                        OGS_FATAL(
                            "Function property: Mismatch of vectorized tensor "
                            "sizes for variable {:s}.",
                            variable_enum_to_string[static_cast<int>(
                                variable)]);
                    }
                }
            };
            std::visit(assign_value,
                       *std::get<VariableArray::DeformationGradient const*>(
                           new_variable_array.address_of(variable)));
        };

        variable_array.visitVariable(
            BaseLib::Overloaded{assign_scalar, assign_kelvin_vector,
                                assign_deformation_gradient},
            variable);
    }
}

static PropertyDataType evaluateExpressions(
    std::vector<Variable> const& variables,
    VariableArray const& new_variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    std::vector<exprtk::expression<double>> const& expressions,
    VariableArray& variable_array, bool const spatial_position_is_required)
{
    std::vector<double> result(expressions.size());

    updateVariableArrayValues(variables, new_variable_array, variable_array);

    std::transform(begin(expressions), end(expressions), begin(result),
                   [t, &pos, spatial_position_is_required](auto const& e)
                   {
                       auto& symbol_table = e.get_symbol_table();
                       symbol_table.get_variable("t")->ref() = t;

                       if (spatial_position_is_required)
                       {
                           if (!pos.getCoordinates())
                           {
                               OGS_FATAL(
                                   "FunctionParameter: The spatial position "
                                   "has to be set by "
                                   "coordinates.");
                           }
                           auto const coords = pos.getCoordinates().value();
                           symbol_table.get_variable("x")->ref() = coords[0];
                           symbol_table.get_variable("y")->ref() = coords[1];
                           symbol_table.get_variable("z")->ref() = coords[2];
                       }

                       return e.value();
                   });

    switch (result.size())
    {
        case 1:
        {
            return result[0];
        }
        case 2:
        {
            return Eigen::Vector2d{result[0], result[1]};
        }
        case 3:
        {
            return Eigen::Vector3d{result[0], result[1], result[2]};
        }
        case 4:
        {
            Eigen::Matrix<double, 2, 2> m;
            m = Eigen::Map<Eigen::Matrix<double, 2, 2> const>(result.data(), 2,
                                                              2);
            return m;
        }
        case 9:
        {
            Eigen::Matrix<double, 3, 3> m;
            m = Eigen::Map<Eigen::Matrix<double, 3, 3> const>(result.data(), 3,
                                                              3);
            return m;
        }
    }
    OGS_FATAL("Cannot convert a vector of size {} to a PropertyDataType",
              result.size());
}

static std::vector<std::string> collectVariables(
    std::vector<std::string> const& value_string_expressions,
    std::vector<std::pair<std::string, std::vector<std::string>>> const&
        dvalue_string_expressions)
{
    std::vector<std::string> variables;

    auto collect_variables = [&](auto string_expressions)
    {
        for (auto const& string_expression : string_expressions)
        {
            if (!exprtk::collect_variables(string_expression, variables))
            {
                OGS_FATAL(
                    "Collecting variables from expression '{}' didn't work.",
                    string_expression);
            }
        }
    };

    collect_variables(value_string_expressions);
    for (auto const& var_name_string_expression : dvalue_string_expressions)
    {
        collect_variables(var_name_string_expression.second);
    }

    BaseLib::makeVectorUnique(variables);
    return variables;
}

Function::Function(
    std::string name,
    std::vector<std::string> const& value_string_expressions,
    std::vector<std::pair<std::string, std::vector<std::string>>> const&
        dvalue_string_expressions,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    name_ = std::move(name);

    auto const variables =
        collectVariables(value_string_expressions, dvalue_string_expressions);

    // filter the strings unrelated to time or spatial position
    static std::unordered_set<std::string> filter_not_variables = {"t", "x",
                                                                   "y", "z"};
    for (auto const& curve : curves)
    {
        filter_not_variables.insert(curve.first);
    }

    variables_ =
        variables |
        std::views::filter([](const std::string& s)
                           { return !filter_not_variables.contains(s); }) |
        std::views::transform([](std::string const& s)
                              { return convertStringToVariable(s); }) |
        ranges::to<std::vector>;

#ifdef _OPENMP
    int const num_threads = omp_get_max_threads();
#else
    int const num_threads = 1;
#endif

    impl2_ = std::make_unique<Implementation<2>>(
        num_threads, variables, value_string_expressions,
        dvalue_string_expressions, curves);
    impl3_ = std::make_unique<Implementation<3>>(
        num_threads, variables, value_string_expressions,
        dvalue_string_expressions, curves);
}

std::variant<Function::Implementation<2>*, Function::Implementation<3>*>
Function::getImplementationForDimensionOfVariableArray(
    VariableArray const& variable_array) const
{
    if (variable_array.is2D())
    {
        return impl2_.get();
    }
    if (variable_array.is3D())
    {
        return impl3_.get();
    }

    OGS_FATAL(
        "Variable array has vectors for 2 and 3 dimensions simultaneously. "
        "Mixed dimensions cannot be dealt within Function evaluation.");
}

PropertyDataType Function::value(VariableArray const& variable_array,
                                 ParameterLib::SpatialPosition const& pos,
                                 double const t, double const /*dt*/) const
{
#ifdef _OPENMP
    int const thread_id = omp_get_thread_num();
#else
    int const thread_id = 0;
#endif
    return std::visit(
        [&](auto&& impl_ptr)
        {
            return evaluateExpressions(variables_, variable_array, pos, t,
                                       impl_ptr->value_expressions[thread_id],
                                       impl_ptr->variable_arrays[thread_id],
                                       impl_ptr->spatial_position_is_required);
        },
        getImplementationForDimensionOfVariableArray(variable_array));
}

PropertyDataType Function::dValue(VariableArray const& variable_array,
                                  Variable const variable,
                                  ParameterLib::SpatialPosition const& pos,
                                  double const t, double const /*dt*/) const
{
#ifdef _OPENMP
    int const thread_id = omp_get_thread_num();
#else
    int const thread_id = 0;
#endif
    return std::visit(
        [&](auto&& impl_ptr)
        {
            auto const it = std::find_if(
                begin(impl_ptr->dvalue_expressions[thread_id]),
                end(impl_ptr->dvalue_expressions[thread_id]),
                [&variable](auto const& v) { return v.first == variable; });

            if (it == end(impl_ptr->dvalue_expressions[thread_id]))
            {
                OGS_FATAL(
                    "Requested derivative with respect to the variable {:s} "
                    "not "
                    "provided for Function-type property {:s}.",
                    variable_enum_to_string[static_cast<int>(variable)], name_);
            }

            return evaluateExpressions(variables_, variable_array, pos, t,
                                       it->second,
                                       impl_ptr->variable_arrays[thread_id],
                                       impl_ptr->spatial_position_is_required);
        },
        getImplementationForDimensionOfVariableArray(variable_array));
}

Function::~Function() = default;

}  // namespace MaterialPropertyLib
