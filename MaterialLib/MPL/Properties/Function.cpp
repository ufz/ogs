// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MaterialLib/MPL/Properties/Function.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>
#include <ranges>
#include <unordered_set>

#include "BaseLib/Algorithm.h"
#include "BaseLib/OgsAsmThreads.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/VectorizedTensor.h"
#include "ParameterLib/ExprtkUtils.h"

namespace MaterialPropertyLib
{

/// Scalar copy operation: transfers one double from a VariableArray into the
/// per-thread symbol table using a precomputed byte offset.
/// The \p variable field is carried only for error messages.
struct ScalarCopyOp
{
    std::ptrdiff_t src_offset;  ///< byte offset from VariableArray start
    double* dst_ptr;            ///< pointer into this thread's symbol table
    Variable variable;          ///< for error messages

    ScalarCopyOp(VariableArray const& base,
                 VariableArray::Scalar* dst,
                 Variable var)
        : src_offset{reinterpret_cast<char const*>(dst) -
                     reinterpret_cast<char const*>(&base)},
          dst_ptr{dst},
          variable{var}
    {
    }

    void apply(VariableArray const& src) const
    {
        double const val = *reinterpret_cast<double const*>(
            reinterpret_cast<char const*>(&src) + src_offset);
        if (std::isnan(val))
        {
            OGS_FATAL(
                "Function property: Scalar variable '{:s}' is not "
                "initialized.",
                variable_enum_to_string[static_cast<int>(variable)]);
        }
        *dst_ptr = val;
    }
};

/// Creates a symbol table for one MPL Function thread, registering standard
/// variables (t, x, y, z), VariableArray fields, and curve wrappers.
/// \p variable_array and \p curve_wrappers are members of the owning
/// PerThreadData and must already be at their final address.
template <int D>
static exprtk::symbol_table<double> createSymbolTable(
    std::vector<std::string> const& expression_symbol_names,
    bool const spatial_position_is_required,
    std::vector<std::string> const& used_curve_names,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    VariableArray& variable_array,
    std::map<std::string, ParameterLib::CurveWrapper>& curve_wrappers)
{
    auto symbol_table =
        ParameterLib::createBaseSymbolTable(spatial_position_is_required);

    std::unordered_set<std::string> curve_name_set(used_curve_names.begin(),
                                                   used_curve_names.end());

    for (auto const& v : expression_symbol_names)
    {
        if (ParameterLib::isBuiltinSymbol(v) || curve_name_set.contains(v))
        {
            continue;
        }

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
                auto& result =
                    address->template emplace<Eigen::Matrix<double, size, 1>>();
                add_vector(result.data(), size);
            },
            [&add_vector](VariableArray::DeformationGradient* address)
            {
                auto constexpr size = MathLib::VectorizedTensor::size(D);
                auto& result =
                    address->template emplace<Eigen::Matrix<double, size, 1>>();
                add_vector(result.data(), size);
            }};

        Variable const variable = convertStringToVariable(v);
        variable_array.visitVariable(add_any_variable, variable);
    }

    ParameterLib::registerCurveWrappers(symbol_table, used_curve_names, curves,
                                        curve_wrappers);
    return symbol_table;
}

/// Symbol table storage and compiled expressions for one OpenMP thread.
/// Each thread receives its own instance so that concurrent evaluations
/// do not race on the exprtk symbol table or the VariableArray scratch space.
/// The symbol table has to point to some fixed addresses. Therefore, it cannot
/// directly use a user-provided VariableArray during evaluation. Instead, a
/// user-provided VariableArray will be copied to this object before expression
/// evaluation.
struct PerThreadData
{
    template <int D>
    PerThreadData(
        std::integral_constant<int, D> /*dim_tag*/,
        std::vector<std::string> const& expression_symbol_names,
        bool spatial_position_is_required,
        std::vector<std::string> const& used_curve_names,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves,
        std::vector<Variable> const& variables_enum,
        std::vector<Variable>* non_scalar_out,
        std::vector<std::string> const& value_string_expressions,
        std::vector<std::pair<std::string, std::vector<std::string>>> const&
            dvalue_string_expressions)
        : symbol_table(createSymbolTable<D>(
              expression_symbol_names, spatial_position_is_required,
              used_curve_names, curves, variable_array, curve_wrappers)),
          symbol_table_cache(symbol_table, spatial_position_is_required)
    {
        buildCopyOps(variables_enum, non_scalar_out);
        compileExpressions(value_string_expressions, dvalue_string_expressions);
    }

    PerThreadData(PerThreadData&& /*other*/)
    {
        OGS_FATAL(
            "MPL's Function property: The internal PerThreadData is not "
            "move-constructible.");
    }
    PerThreadData& operator=(PerThreadData&& other) noexcept = delete;
    PerThreadData(PerThreadData const&) = delete;
    PerThreadData& operator=(PerThreadData const&) = delete;

    /// Curve wrappers owned by this thread; must outlive the symbol table.
    std::map<std::string, ParameterLib::CurveWrapper> curve_wrappers;

    /// Storage for vectorial variables (KelvinVector / DeformationGradient).
    /// Must be declared before symbol_table so it is initialized first
    /// (symbol_table stores pointers into variable_array).
    VariableArray variable_array;

    /// Symbol table holding variable storage. Must be destroyed after
    /// expressions (exprtk requirement) and after curve_wrappers and
    /// variable_array.
    exprtk::symbol_table<double> symbol_table;

    /// Value expressions evaluating the property value.
    /// Multi-component expressions corresponding one-to-one with the input
    /// XML entries, e.g.:
    /// ```xml
    /// <value>
    ///   <expression>f1</expression>
    ///   <expression>f2</expression>
    /// </value>
    /// ```
    std::vector<exprtk::expression<double>> value_expressions;

    /// Derivative expressions w.r.t. a differentiation variable.
    /// Outer vector is an unordered set of (Variable, expressions) pairs.
    /// Each pair's second element holds multi-component expressions,
    /// corresponding one-to-one with the input XML entries, e.g.:
    /// ```xml
    /// <dvalue>
    ///   <variable>T</variable>
    ///   <expression>df1/dT</expression>
    ///   <expression>df2/dT</expression>
    /// </dvalue>
    /// ```
    std::vector<std::pair<Variable, std::vector<exprtk::expression<double>>>>
        dvalue_expressions;

    /// Cached pointers to symbol table variables (t, x, y, z).
    ParameterLib::SymbolTableCache symbol_table_cache;

    /// Scalar copy operations into this thread's symbol table.
    std::vector<ScalarCopyOp> scalar_copy_ops;

    /// Update non-scalar (KelvinVector / DeformationGradient) variables in
    /// this thread's variable_array from \p new_variable_array.
    void updateNonScalarVariables(
        std::vector<Variable> const& non_scalar_variables,
        VariableArray const& new_variable_array);

    /// Evaluate \p expressions against the given variable array and position,
    /// returning the result as a PropertyDataType.
    /// \param non_scalar_variables  Non-scalar variables from Implementation
    ///                              (shared read-only).
    /// \param new_variable_array    Input variable array for this evaluation.
    /// \param pos                   Spatial position (coordinates required when
    ///                              spatial variables x, y, z are used).
    /// \param t                     Current time.
    /// \param expressions           Compiled exprtk expressions to evaluate.
    PropertyDataType evaluate(
        std::vector<Variable> const& non_scalar_variables,
        VariableArray const& new_variable_array,
        ParameterLib::SpatialPosition const& pos, double t,
        std::vector<exprtk::expression<double>> const& expressions);

    /// Build scalar_copy_ops for this thread.
    /// \param variables_enum    All variables to visit.
    /// \param non_scalar_out    If non-null, non-scalar variables are appended
    ///                          here (pass only for thread 0).
    void buildCopyOps(std::vector<Variable> const& variables_enum,
                      std::vector<Variable>* non_scalar_out)
    {
        for (auto const variable : variables_enum)
        {
            auto const add_non_scalar = [&]()
            {
                if (non_scalar_out)
                {
                    non_scalar_out->push_back(variable);
                }
            };

            variable_array.visitVariable(
                BaseLib::Overloaded{[&](VariableArray::Scalar* dst)
                                    {
                                        scalar_copy_ops.emplace_back(
                                            variable_array, dst, variable);
                                    },
                                    [&](VariableArray::KelvinVector*)
                                    { add_non_scalar(); },
                                    [&](VariableArray::DeformationGradient*)
                                    { add_non_scalar(); }},
                variable);
        }
    }

    /// Compile value and dValue expressions for this thread.
    /// \param value_string_expressions  Raw value expression strings.
    /// \param dvalue_string_expressions Variable name paired with its
    ///                                  derivative expression strings.
    void compileExpressions(
        std::vector<std::string> const& value_string_expressions,
        std::vector<std::pair<std::string, std::vector<std::string>>> const&
            dvalue_string_expressions)
    {
        value_expressions = ParameterLib::compileExpressions(
            symbol_table, value_string_expressions);

        for (auto const& [variable_name, string_expressions] :
             dvalue_string_expressions)
        {
            if (string_expressions.size() != value_string_expressions.size())
            {
                OGS_FATAL(
                    "The number of dValue expressions ({:d}) for variable "
                    "'{:s}' does not match the number of value expressions "
                    "({:d}).",
                    string_expressions.size(), variable_name,
                    value_string_expressions.size());
            }
            dvalue_expressions.emplace_back(
                convertStringToVariable(variable_name),
                ParameterLib::compileExpressions(symbol_table,
                                                 string_expressions));
        }
    }
};

template <int D>
struct Function::Implementation
{
    using Expression = exprtk::expression<double>;

    Implementation(
        int num_threads,
        std::vector<std::string> const& expression_symbol_names,
        std::vector<Variable> const& variables_enum,
        std::vector<std::string> const& value_string_expressions,
        std::vector<std::pair<std::string, std::vector<std::string>>> const&
            dvalue_string_expressions,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves);

    /// Per-thread data; indexed by omp_get_thread_num().
    std::vector<PerThreadData> per_thread_data;

    /// Non-scalar variables (KelvinVector / DeformationGradient); shared
    /// read-only across all threads.
    std::vector<Variable> non_scalar_variables;

    bool spatial_position_is_required = false;
};

template <int D>
Function::Implementation<D>::Implementation(
    int num_threads,
    std::vector<std::string> const& expression_symbol_names,
    std::vector<Variable> const& variables_enum,
    std::vector<std::string> const& value_string_expressions,
    std::vector<std::pair<std::string, std::vector<std::string>>> const&
        dvalue_string_expressions,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    spatial_position_is_required =
        ParameterLib::hasSpatialVariables(expression_symbol_names);

    // Determine which curves are actually referenced in the expressions.
    // The expression symbol names list (from exprtk::collect_variables)
    // contains both VariableArray field names and curve names as plain
    // identifiers.
    auto const used_curve_names =
        ParameterLib::collectUsedCurveNames(expression_symbol_names, curves);

    per_thread_data.reserve(num_threads);
    for (int thread_id = 0; thread_id < num_threads; ++thread_id)
    {
        per_thread_data.emplace_back(
            std::integral_constant<int, D>{}, expression_symbol_names,
            spatial_position_is_required, used_curve_names, curves,
            variables_enum, thread_id == 0 ? &non_scalar_variables : nullptr,
            value_string_expressions, dvalue_string_expressions);
    }
}

void PerThreadData::updateNonScalarVariables(
    std::vector<Variable> const& non_scalar_variables,
    VariableArray const& new_variable_array)
{
    for (auto const& variable : non_scalar_variables)
    {
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
            BaseLib::Overloaded{
                [](VariableArray::Scalar*)
                {
                    OGS_FATAL(
                        "Function property: updateNonScalarVariables called "
                        "with a scalar variable.");
                },
                assign_kelvin_vector, assign_deformation_gradient},
            variable);
    }
}

/// Helper to get fixed-size array from expressions for known sizes.
/// This avoids heap allocation.
template <std::size_t N>
static std::array<double, N> evaluateToArray(
    std::vector<exprtk::expression<double>> const& expressions)
{
    std::array<double, N> result{};
    for (std::size_t i = 0; i < N; ++i)
    {
        result[i] = expressions[i].value();
    }
    return result;
}

PropertyDataType PerThreadData::evaluate(
    std::vector<Variable> const& non_scalar_variables,
    VariableArray const& new_variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    std::vector<exprtk::expression<double>> const& expressions)
{
    // Fast path: direct offset arithmetic.
    for (auto const& op : scalar_copy_ops)
    {
        op.apply(new_variable_array);
    }
    // Fallback for KelvinVector / DeformationGradient variables.
    updateNonScalarVariables(non_scalar_variables, new_variable_array);

    // Set symbol table variables using cached pointers.
    symbol_table_cache.setTimeAndPosition(t, pos);

    // Dispatch on expression count: template specialisations fix N at compile
    // time, enabling stack allocation in evaluateToArray and branch-free
    // conversion in fromArray.
    auto const n = expressions.size();
    switch (n)
    {
        case 1:
            return fromArray(evaluateToArray<1>(expressions));
        case 2:
            return fromArray(evaluateToArray<2>(expressions));
        case 3:
            return fromArray(evaluateToArray<3>(expressions));
        case 4:
            return fromArray(evaluateToArray<4>(expressions));
        case 6:
            return fromArray(evaluateToArray<6>(expressions));
        case 9:
            return fromArray(evaluateToArray<9>(expressions));
        default:
            OGS_FATAL(
                "Cannot convert a vector of size {} to a PropertyDataType", n);
    }
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

    // Collect variables from all expressions (value and dvalue) so the symbol
    // table is complete even when a variable appears only in a derivative.
    auto all_exprs = value_string_expressions;
    for (auto const& [_, exprs] : dvalue_string_expressions)
    {
        all_exprs.insert(all_exprs.end(), exprs.begin(), exprs.end());
    }
    auto const expression_symbol_names =
        ParameterLib::collectVariables(all_exprs);

    required_variables_enum_ =
        expression_symbol_names |
        std::views::filter(
            [&curves](std::string const& s)
            {
                return !ParameterLib::isBuiltinSymbol(s) && !curves.contains(s);
            }) |
        std::views::transform([](std::string const& s)
                              { return convertStringToVariable(s); }) |
        ranges::to<std::vector>;

    int const num_threads = BaseLib::getNumberOfThreads();

    impl2_ = std::make_unique<Implementation<2>>(
        num_threads, expression_symbol_names, required_variables_enum_,
        value_string_expressions, dvalue_string_expressions, curves);
    impl3_ = std::make_unique<Implementation<3>>(
        num_threads, expression_symbol_names, required_variables_enum_,
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
            if (thread_id >= static_cast<int>(impl_ptr->per_thread_data.size()))
            {
                OGS_FATAL(
                    "In Function-type property '{:s}' evaluation the "
                    "OMP-thread with id {:d} exceeds the number of allocated "
                    "threads {:d}.",
                    name_, thread_id, impl_ptr->per_thread_data.size());
            }
            auto& thread_data = impl_ptr->per_thread_data[thread_id];
            return thread_data.evaluate(impl_ptr->non_scalar_variables,
                                        variable_array, pos, t,
                                        thread_data.value_expressions);
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
            if (thread_id >= static_cast<int>(impl_ptr->per_thread_data.size()))
            {
                OGS_FATAL(
                    "In Function-type property '{:s}' evaluation the "
                    "OMP-thread with id {:d} exceeds the number of allocated "
                    "threads {:d}.",
                    name_, thread_id, impl_ptr->per_thread_data.size());
            }
            auto& thread_data = impl_ptr->per_thread_data[thread_id];
            auto const it = std::ranges::find_if(
                thread_data.dvalue_expressions,
                [&variable](auto const& v) { return v.first == variable; });

            if (it == end(thread_data.dvalue_expressions))
            {
                OGS_FATAL(
                    "Requested derivative with respect to the variable {:s} "
                    "not provided for Function-type property {:s}.",
                    variable_enum_to_string[static_cast<int>(variable)], name_);
            }

            return thread_data.evaluate(impl_ptr->non_scalar_variables,
                                        variable_array, pos, t, it->second);
        },
        getImplementationForDimensionOfVariableArray(variable_array));
}

Function::~Function() = default;

}  // namespace MaterialPropertyLib
