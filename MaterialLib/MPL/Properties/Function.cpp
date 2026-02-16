// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
#include "BaseLib/OgsAsmThreads.h"
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

struct SymbolTableCache
{
    double* t_ptr = nullptr;
    double* x_ptr = nullptr;
    double* y_ptr = nullptr;
    double* z_ptr = nullptr;
};

struct ScalarCopyOp
{
    std::ptrdiff_t src_offset;  ///< byte offset from VariableArray start
    double* dst_ptr;            ///< pointer into this thread's variable_array
    Variable variable;          ///< for error messages
};

template <int D>
struct Function::Implementation
{
    using Expression = exprtk::expression<double>;

    Implementation(
        int num_threads,
        std::vector<std::string> const& variables,
        std::vector<Variable> const& variables_enum,
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

    /// Create symbol tables for all threads and populate scalar_copy_ops.
    std::vector<exprtk::symbol_table<double>> createSymbolTables(
        int num_threads,
        std::vector<std::string> const& variables,
        std::vector<Variable> const& variables_enum,
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

    /// Cached pointers to symbol table variables (t, x, y, z) per thread.
    /// Eliminates string-based lookups in get_variable().
    std::vector<SymbolTableCache> symbol_table_caches;

    /// Per-thread list of scalar copy operations (src offset → dst pointer).
    /// Built once at construction; used instead of updateVariableArrayValues.
    std::vector<std::vector<ScalarCopyOp>> scalar_copy_ops;

    /// Non-scalar variables (KelvinVector / DeformationGradient) that still
    /// go through the generic updateVariableArrayValues path.
    std::vector<Variable> non_scalar_variables;

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
        if (v == "x" || v == "y" || v == "z")
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
    std::vector<Variable> const& variables_enum,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    std::vector<exprtk::symbol_table<double>> symbol_tables(num_threads);
    symbol_table_caches.resize(num_threads);

    for (int thread_id = 0; thread_id < num_threads; ++thread_id)
    {
        symbol_tables[thread_id] =
            createSymbolTable(variables, curves, variable_arrays[thread_id]);

        // Cache pointers to frequently accessed variables to avoid
        // std::map lookups in get_variable().
        symbol_table_caches[thread_id].t_ptr =
            &(symbol_tables[thread_id].get_variable("t")->ref());

        if (spatial_position_is_required)
        {
            symbol_table_caches[thread_id].x_ptr =
                &(symbol_tables[thread_id].get_variable("x")->ref());
            symbol_table_caches[thread_id].y_ptr =
                &(symbol_tables[thread_id].get_variable("y")->ref());
            symbol_table_caches[thread_id].z_ptr =
                &(symbol_tables[thread_id].get_variable("z")->ref());
        }
    }

    // Build per-thread scalar copy ops (and shared non-scalar fallback list).
    scalar_copy_ops.resize(num_threads);
    for (int thread_id = 0; thread_id < num_threads; ++thread_id)
    {
        for (auto const variable : variables_enum)
        {
            variable_arrays[thread_id].visitVariable(
                BaseLib::Overloaded{
                    [&](VariableArray::Scalar* dst)
                    {
                        auto const src_offset =
                            reinterpret_cast<char const*>(dst) -
                            reinterpret_cast<char const*>(
                                &variable_arrays[thread_id]);
                        scalar_copy_ops[thread_id].push_back(
                            {src_offset, dst, variable});
                    },
                    [&](VariableArray::KelvinVector*)
                    {
                        if (thread_id == 0)
                        {
                            non_scalar_variables.push_back(variable);
                        }
                    },
                    [&](VariableArray::DeformationGradient*)
                    {
                        if (thread_id == 0)
                        {
                            non_scalar_variables.push_back(variable);
                        }
                    }},
                variable);
        }
    }

    return symbol_tables;
}

template <int D>
Function::Implementation<D>::Implementation(
    int num_threads,
    std::vector<std::string> const& variables,
    std::vector<Variable> const& variables_enum,
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

    auto symbol_tables =
        createSymbolTables(num_threads, variables, variables_enum, curves);

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

template <std::size_t N>
static PropertyDataType evaluateExpressionsImpl(
    std::vector<ScalarCopyOp> const& scalar_copy_ops,
    std::vector<Variable> const& non_scalar_variables,
    VariableArray const& new_variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    std::vector<exprtk::expression<double>> const& expressions,
    VariableArray& variable_array, bool const spatial_position_is_required,
    SymbolTableCache const& cache)
{
    // Fast path: direct offset arithmetic, no switches, no variant visits.
    for (auto const& op : scalar_copy_ops)
    {
        double const val = *reinterpret_cast<double const*>(
            reinterpret_cast<char const*>(&new_variable_array) + op.src_offset);
        if (std::isnan(val))
        {
            OGS_FATAL(
                "Function property: Scalar variable '{:s}' is not "
                "initialized.",
                variable_enum_to_string[static_cast<int>(op.variable)]);
        }
        *op.dst_ptr = val;
    }
    // Fallback for KelvinVector / DeformationGradient variables.
    if (!non_scalar_variables.empty())
    {
        updateVariableArrayValues(non_scalar_variables, new_variable_array,
                                  variable_array);
    }

    // Use std::array for stack allocation (no heap allocation).
    std::array<double, N> result{};

    // Set symbol table variables using cached pointers (eliminates
    // std::map lookups).
    *cache.t_ptr = t;

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
        *cache.x_ptr = coords[0];
        *cache.y_ptr = coords[1];
        *cache.z_ptr = coords[2];
    }

    for (std::size_t i = 0; i < N; ++i)
    {
        result[i] = expressions[i].value();
    }

    // Convert result to appropriate PropertyDataType based on size.
    if constexpr (N == 1)
    {
        return result[0];
    }
    else if constexpr (N == 2)
    {
        return Eigen::Vector2d{result[0], result[1]};
    }
    else if constexpr (N == 3)
    {
        return Eigen::Vector3d{result[0], result[1], result[2]};
    }
    else if constexpr (N == 4)
    {
        Eigen::Matrix<double, 2, 2> m;
        m = Eigen::Map<Eigen::Matrix<double, 2, 2> const>(result.data(), 2, 2);
        return m;
    }
    else if constexpr (N == 6)
    {
        OGS_FATAL("Cannot convert a vector of size {} to a PropertyDataType",
                  N);
    }
    else if constexpr (N == 9)
    {
        Eigen::Matrix<double, 3, 3> m;
        m = Eigen::Map<Eigen::Matrix<double, 3, 3> const>(result.data(), 3, 3);
        return m;
    }
}

static PropertyDataType evaluateExpressions(
    std::vector<ScalarCopyOp> const& scalar_copy_ops,
    std::vector<Variable> const& non_scalar_variables,
    VariableArray const& new_variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    std::vector<exprtk::expression<double>> const& expressions,
    VariableArray& variable_array, bool const spatial_position_is_required,
    SymbolTableCache const& cache)
{
    switch (expressions.size())
    {
        case 1:
            return evaluateExpressionsImpl<1>(
                scalar_copy_ops, non_scalar_variables, new_variable_array, pos,
                t, expressions, variable_array, spatial_position_is_required,
                cache);
        case 2:
            return evaluateExpressionsImpl<2>(
                scalar_copy_ops, non_scalar_variables, new_variable_array, pos,
                t, expressions, variable_array, spatial_position_is_required,
                cache);
        case 3:
            return evaluateExpressionsImpl<3>(
                scalar_copy_ops, non_scalar_variables, new_variable_array, pos,
                t, expressions, variable_array, spatial_position_is_required,
                cache);
        case 4:
            return evaluateExpressionsImpl<4>(
                scalar_copy_ops, non_scalar_variables, new_variable_array, pos,
                t, expressions, variable_array, spatial_position_is_required,
                cache);
        case 6:
            return evaluateExpressionsImpl<6>(
                scalar_copy_ops, non_scalar_variables, new_variable_array, pos,
                t, expressions, variable_array, spatial_position_is_required,
                cache);
        case 9:
            return evaluateExpressionsImpl<9>(
                scalar_copy_ops, non_scalar_variables, new_variable_array, pos,
                t, expressions, variable_array, spatial_position_is_required,
                cache);
        default:
            OGS_FATAL(
                "Cannot convert a vector of size {} to a PropertyDataType",
                expressions.size());
    }
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

    required_variables_enum_ =
        variables |
        std::views::filter([](const std::string& s)
                           { return !filter_not_variables.contains(s); }) |
        std::views::transform([](std::string const& s)
                              { return convertStringToVariable(s); }) |
        ranges::to<std::vector>;

    auto const get_number_omp_threads = []()
    {
#ifdef _OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    };

    int const num_threads = std::max(BaseLib::getNumberOfAssemblyThreads(),
                                     get_number_omp_threads());

    impl2_ = std::make_unique<Implementation<2>>(
        num_threads, variables, required_variables_enum_,
        value_string_expressions, dvalue_string_expressions, curves);
    impl3_ = std::make_unique<Implementation<3>>(
        num_threads, variables, required_variables_enum_,
        value_string_expressions, dvalue_string_expressions, curves);
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
            if (thread_id >=
                static_cast<int>(impl_ptr->value_expressions.size()))
            {
                OGS_FATAL(
                    "In Function-type property '{:s}' evaluation the "
                    "OMP-thread with id {:d} exceeds the number of allocated "
                    "threads {:d}.",
                    name_, thread_id, impl_ptr->value_expressions.size());
            }
            return evaluateExpressions(
                impl_ptr->scalar_copy_ops[thread_id],
                impl_ptr->non_scalar_variables, variable_array, pos, t,
                impl_ptr->value_expressions[thread_id],
                impl_ptr->variable_arrays[thread_id],
                impl_ptr->spatial_position_is_required,
                impl_ptr->symbol_table_caches[thread_id]);
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
            if (thread_id >=
                static_cast<int>(impl_ptr->dvalue_expressions.size()))
            {
                OGS_FATAL(
                    "In Function-type property '{:s}' evaluation the "
                    "OMP-thread with id {:d} exceeds the number of allocated "
                    "threads {:d}.",
                    name_, thread_id, impl_ptr->value_expressions.size());
            }
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

            return evaluateExpressions(
                impl_ptr->scalar_copy_ops[thread_id],
                impl_ptr->non_scalar_variables, variable_array, pos, t,
                it->second, impl_ptr->variable_arrays[thread_id],
                impl_ptr->spatial_position_is_required,
                impl_ptr->symbol_table_caches[thread_id]);
        },
        getImplementationForDimensionOfVariableArray(variable_array));
}

Function::~Function() = default;

}  // namespace MaterialPropertyLib
