/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MaterialLib/MPL/Properties/Function.h"

#include <exprtk.hpp>
#include <numeric>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>

#include "BaseLib/Algorithm.h"

namespace MaterialPropertyLib
{
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

class Function::Implementation
{
public:
    using Expression = exprtk::expression<double>;

public:
    Implementation(
        std::vector<std::string> const& variables,
        std::vector<std::string> const& value_string_expressions,
        std::vector<std::pair<std::string, std::vector<std::string>>> const&
            dvalue_string_expressions);

private:
    /// Create symbol table for given variables and populates the variable_array
    /// as needed.
    exprtk::symbol_table<double> createSymbolTable(
        std::vector<std::string> const& variables);

public:
    /// Value expressions.
    /// Multiple expressions are representing vector-valued functions.
    std::vector<Expression> value_expressions;

    /// Derivative expressions with respect to the variable.
    /// Multiple expressions are representing vector-valued functions.
    std::vector<std::pair<Variable, std::vector<Expression>>>
        dvalue_expressions;

    /// Stores values for evaluation of vectorial quantities. Needed for
    /// constant pointers for exprtk.
    mutable VariableArray variable_array;
};

exprtk::symbol_table<double> Function::Implementation::createSymbolTable(
    std::vector<std::string> const& variables)
{
    exprtk::symbol_table<double> symbol_table;

    for (auto const& v : variables)
    {
        auto add_scalar = [&v, &symbol_table](double& value)
        { symbol_table.add_variable(v, value); };

        Variable const variable = convertStringToVariable(v);

        auto add_any_variable = [&add_scalar, &variable]<typename T>(T* address)
        {
            if constexpr (std::is_same_v<VariableArray::Scalar, T>)
            {
                add_scalar(*address);
            }
            else
            {
                OGS_FATAL(
                    "Function property: Non-scalar quantities (for variable "
                    "{:s}) are not handled by the implementation yet.",
                    variable_enum_to_string[static_cast<int>(variable)]);
            }
        };

        std::visit(add_any_variable, variable_array.address_of(variable));
    }
    return symbol_table;
}

Function::Implementation::Implementation(
    std::vector<std::string> const& variables,
    std::vector<std::string> const& value_string_expressions,
    std::vector<std::pair<std::string, std::vector<std::string>>> const&
        dvalue_string_expressions)
{
    auto symbol_table = createSymbolTable(variables);

    // value expressions.
    value_expressions =
        compileExpressions(symbol_table, value_string_expressions);

    // dValue expressions.
    for (auto const& [variable_name, string_expressions] :
         dvalue_string_expressions)
    {
        dvalue_expressions.emplace_back(
            convertStringToVariable(variable_name),
            compileExpressions(symbol_table, string_expressions));
    }
}

static void updateVariableArrayValues(std::vector<Variable> const& variables,
                                      VariableArray const& new_variable_array,
                                      VariableArray& variable_array)
{
    for (auto const& variable : variables)
    {
        auto assign_variable =
            [&variable, &new_variable_array]<typename T>(T* address)
        {
            if constexpr (std::is_same_v<VariableArray::Scalar, T>)
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
            }
            else
            {
                OGS_FATAL(
                    "Function property: Non-scalar quantities (for variable "
                    "{:s}) are not handled by the implementation yet.",
                    variable_enum_to_string[static_cast<int>(variable)]);
            }
        };
        std::visit(assign_variable, variable_array.address_of(variable));
    }
}

static PropertyDataType evaluateExpressions(
    std::vector<Variable> const& variables,
    VariableArray const& new_variable_array,
    std::vector<exprtk::expression<double>> const& expressions,
    VariableArray& variable_array,
    std::mutex& mutex)
{
    std::vector<double> result(expressions.size());

    {
        std::lock_guard lock_guard(mutex);
        updateVariableArrayValues(variables, new_variable_array,
                                  variable_array);

        std::transform(begin(expressions), end(expressions), begin(result),
                       [](auto const& e) { return e.value(); });
    }

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
        dvalue_string_expressions)
{
    name_ = std::move(name);

    auto const variables =
        collectVariables(value_string_expressions, dvalue_string_expressions);
    variables_ =
        variables |
        ranges::views::transform([](std::string const& s)
                                 { return convertStringToVariable(s); }) |
        ranges::to<std::vector>;

    impl_ptr_ = std::make_unique<Implementation>(
        variables, value_string_expressions, dvalue_string_expressions);
}

PropertyDataType Function::value(VariableArray const& variable_array,
                                 ParameterLib::SpatialPosition const& /*pos*/,
                                 double const /*t*/, double const /*dt*/) const
{
    return evaluateExpressions(variables_, variable_array,
                               impl_ptr_->value_expressions,
                               impl_ptr_->variable_array, mutex_);
}

PropertyDataType Function::dValue(VariableArray const& variable_array,
                                  Variable const variable,
                                  ParameterLib::SpatialPosition const& /*pos*/,
                                  double const /*t*/, double const /*dt*/) const
{
    auto const it = std::find_if(begin(impl_ptr_->dvalue_expressions),
                                 end(impl_ptr_->dvalue_expressions),
                                 [&variable](auto const& v)
                                 { return v.first == variable; });

    if (it == end(impl_ptr_->dvalue_expressions))
    {
        OGS_FATAL(
            "Requested derivative with respect to the variable {:s} not "
            "provided for Function-type property {:s}.",
            variable_enum_to_string[static_cast<int>(variable)], name_);
    }

    return evaluateExpressions(variables_, variable_array, it->second,
                               impl_ptr_->variable_array, mutex_);
}

Function::~Function() = default;

}  // namespace MaterialPropertyLib
