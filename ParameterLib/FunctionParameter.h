/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <utility>
#include <vector>

#include <exprtk.hpp>

#include "Parameter.h"
#include "Utils.h"

namespace ParameterLib
{
/// A parameter class evaluating functions defined by
/// user-provided mathematical expressions.
///
/// Currently, x, y, and z are supported as variables
/// of the functions.
template <typename T>
struct FunctionParameter final : public Parameter<T>
{
    using symbol_table_t = exprtk::symbol_table<T>;
    using expression_t = exprtk::expression<T>;
    using parser_t = exprtk::parser<T>;
    using error_t = exprtk::parser_error::type;

    /**
     * Constructing from a vector of expressions
     *
     * @param name        the parameter's name
     * @param vec_expression_str  a vector of mathematical expressions
     * The vector size specifies the number of components of the parameter.
     */
    FunctionParameter(std::string const& name,
                      std::vector<std::string> const& vec_expression_str)
        : Parameter<T>(name, nullptr), _vec_expression_str(vec_expression_str)
    {
        _symbol_table.add_constants();
        _symbol_table.create_variable("x");
        _symbol_table.create_variable("y");
        _symbol_table.create_variable("z");

        _vec_expression.resize(_vec_expression_str.size());
        for (unsigned i = 0; i < _vec_expression_str.size(); i++)
        {
            _vec_expression[i].register_symbol_table(_symbol_table);
            parser_t parser;
            if (!parser.compile(_vec_expression_str[i], _vec_expression[i]))
            {
                OGS_FATAL("Error: {:s}\tExpression: {:s}\n",
                          parser.error(),
                          _vec_expression_str[i]);
            }
        }
    }

    bool isTimeDependent() const override { return false; }

    int getNumberOfComponents() const override
    {
        return _vec_expression.size();
    }

    std::vector<T> operator()(double const /*t*/,
                              SpatialPosition const& pos) const override
    {
        std::vector<T> cache(getNumberOfComponents());
        auto& x = _symbol_table.get_variable("x")->ref();
        auto& y = _symbol_table.get_variable("y")->ref();
        auto& z = _symbol_table.get_variable("z")->ref();
        if (!pos.getCoordinates())
        {
            OGS_FATAL(
                "FunctionParameter: The spatial position has to be set by "
                "coordinates.");
        }
        auto const coords = pos.getCoordinates().get();
        x = coords[0];
        y = coords[1];
        z = coords[2];

        for (unsigned i = 0; i < _vec_expression.size(); i++)
        {
            cache[i] = _vec_expression[i].value();
        }

        if (!this->_coordinate_system)
        {
            return cache;
        }

        return this->rotateWithCoordinateSystem(cache, pos);
    }

private:
    std::vector<std::string> const _vec_expression_str;
    symbol_table_t _symbol_table;
    std::vector<expression_t> _vec_expression;
};

std::unique_ptr<ParameterBase> createFunctionParameter(
    std::string const& name, BaseLib::ConfigTree const& config);

}  // namespace ParameterLib
