/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MaterialLib/MPL/Properties/Function.h"

#include <numeric>

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
    expressions.resize(string_expressions.size());
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

static void updateVariableValues(
    std::vector<std::pair<int, double*>> const& symbol_values,
    VariableArray const& variable_array)
{
    for (auto& index_value_ptr_pair : symbol_values)
    {
        auto const index = index_value_ptr_pair.first;

        double* value_ptr = index_value_ptr_pair.second;
        std::visit(
            [&value_ptr, &index](auto&& v)
            {
                using T = std::decay_t<decltype(v)>;
                if constexpr (std::is_same_v<T, std::monostate>)
                {
                    OGS_FATAL(
                        "Function property: variable {:s} value needed for "
                        "evaluation of the expression was not set by the "
                        "caller.",
                        variable_enum_to_string[index]);
                }
                else if constexpr (std::is_same_v<T, double>)
                {
                    *value_ptr = v;
                }
                else
                {
                    OGS_FATAL(
                        "Function property: not implemented handling for a "
                        "type {:s} of variable {:s}.",
                        typeid(T).name(), variable_enum_to_string[index]);
                }
            },
            variable_array[index]);
    }
}

static PropertyDataType evaluateExpressions(
    std::vector<std::pair<int, double*>> const& symbol_values,
    VariableArray const& variable_array,
    std::vector<exprtk::expression<double>> const& expressions)
{
    updateVariableValues(symbol_values, variable_array);

    std::vector<double> result(expressions.size());
    std::transform(begin(expressions), end(expressions), begin(result),
                   [](auto const& e) { return e.value(); });

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
                OGS_FATAL("Collecting variables didn't work.");
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

    // Create symbol table for used variables.
    exprtk::symbol_table<double> symbol_table;

    for (auto const& v :
         collectVariables(value_string_expressions, dvalue_string_expressions))
    {
        symbol_table.create_variable(v);
        // Store variables index in the variable array and the pointer to the
        // value in the symbol table for fast access later.
        int const variable_array_index =
            static_cast<int>(convertStringToVariable(v));
        symbol_values_.emplace_back(variable_array_index,
                                    &symbol_table.get_variable(v)->ref());
    }

    // value expressions.
    value_expressions_ =
        compileExpressions(symbol_table, value_string_expressions);

    // dValue expressions.
    for (auto const& [variable_name, string_expressions] :
         dvalue_string_expressions)
    {
        dvalue_expressions_.emplace_back(
            convertStringToVariable(variable_name),
            compileExpressions(symbol_table, string_expressions));
    }
}

PropertyDataType Function::value(VariableArray const& variable_array,
                                 ParameterLib::SpatialPosition const& /*pos*/,
                                 double const /*t*/, double const /*dt*/) const
{
    return evaluateExpressions(symbol_values_, variable_array,
                               value_expressions_);
}

PropertyDataType Function::dValue(VariableArray const& variable_array,
                                  Variable const primary_variable,
                                  ParameterLib::SpatialPosition const& /*pos*/,
                                  double const /*t*/, double const /*dt*/) const
{
    auto const it = std::find_if(begin(dvalue_expressions_),
                                 end(dvalue_expressions_),
                                 [&primary_variable](auto const& v)
                                 { return v.first == primary_variable; });

    if (it == end(dvalue_expressions_))
    {
        OGS_FATAL(
            "Requested derivative with respect to the variable {:s} not "
            "provided for Function-type property {:s}.",
            variable_enum_to_string[static_cast<int>(primary_variable)], name_);
    }

    return evaluateExpressions(symbol_values_, variable_array, it->second);
}

}  // namespace MaterialPropertyLib
